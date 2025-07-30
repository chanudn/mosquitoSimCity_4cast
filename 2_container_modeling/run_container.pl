#!/usr/bin/perl
use warnings;
use strict;
use Statistics::Descriptive;
use File::Spec; # for concatenating filepaths safely

#File: run_container.pl
#Description: Runs container modeling program for specified experiments.
#             Constructs parameter file based on a template, then runs program.
#Input: Container program, specify parameters here, input met data
#Output: Container program output files
#Dan Steinhoff
# 13 March 2013

# Revision history (CNY)
#
# 2023-??-??, CNY: Added more temperature bins for frequency stats.
#                  Added a run name that summarizes a few specs: city, timespan, manual fill
#                  (this run name is also used for generating directory name).
#                  Set whs_mul (initial water height multiplier) to 0.1 (originally 0.0).
#                  Note: Set cflag (whether clouds specified) to "F" when using GEOS S2S input.
#                  Note: Set sflag (whether soil temp specified) to "F" when using GEOS S2S input.
# 2023-08-09, CNY: Added manual emptying.
# 2023-08-31, CNY: Major rewrite.
#                    - Removed plotting.
#                    - Created helper function lookup_vars() so that a small number of basic user specifications
#                      can be used to create numerous other specifications. This simplifies specification of
#                      variables that are linked to one another (e.g., lat/lon values are looked up by the helper
#                      fn based on the user-specified city name, rather than lat/lon being user-specified as well).
#                    - Consolidated the main body of code into the new function whatchem_main(),
#                      making it easier to run the model over a suite of user specifications.
# 2023-09-05, CNY: Removed the initial water and ground temps (tw, tg) as user specifications. These are now 
#                  estimated as the initial air temperature from the input climate file. Note: tg is ignored if 
#                  the input climate file includes soil temperature data (i.e., if sflag is TRUE).
# 2023-09-06, CNY: Fixed bug: added sprintf() formatting for $tw, $tg, and $tStep_model to match the original code.
#                  Without this formatting the Fortran code interprets these values incorrectly.
# 2023-10-28, CNY: Wrote main code for running model (FOR loops that span the parameter space of interest).
# 2025-05-14, CNY: Added handling of GEOS input data (and updated filepaths to point to whatchem-2.1_CNY_4cast/).
#                  Simplified model-running loop via helper fns is_leap_year(), get_days_in_month(), make_timespan().
#                  Deleted code for climatology model runs since it's unused.


#### USER SPECIFICATIONS ####################################################################

# Container parameters
# These correspond to more specific definitions found in the lookup_vars() function
# (e.g., the $name_refl value here corresponds to a specific reflectivity defined in lookup_vars()).
my $name_cont  = "Bucket";           # Name of container
my $name_refl  = "Gray";             # Name of reflectivity (color)
my $name_shade = "HalfShade";        # Name of shade condition

# Water conditions
my $wh_mul = 0.1;     # Multiplier of initial (water height / container height) ratio

# Flags
my $doClim = "F";      # Whether input data corresponds to climatology (T/F) (ADDED ON 2023-09-10, CNY)
my $cflag  = "T";      # Whether cloud fraction data is input (T/F)
my $sflag  = "T";      # Whether soil temperature data is input (T/F)
my $mfill  = "F";      # Manual fill flag (T/F)
my $fill_intv = 168;   # Manual fill interval, in hours (must be integer)
my $mempty = "F";      # Manual empty flag (T/F)                           (ADDED ON 2023-08-09, CNY)
my $empty_intv = 168;  # Manual empty interval, in hours (must be integer) (ADDED ON 2023-08-09, CNY)

# Misc.
my $tStep_model = 60.00;    # Model time step in seconds
my $tStep_input = 3600.00;  # Time interval of input data in seconds
my $missVal = -9999.00;     # Missing value of input data

# Filepaths
my $path_exe        = "/home/local/WIN/cyasana1/Code/project_whatchem/whatchem-2.1_CNY_4cast/2_container_modeling/src"; # Path to container.exe file
my $dir_input_base  = "/mnt/redwood/local_drive/chanud/whatchem_2.1_CNY_4cast/input_data";                              # Directory of input data
my $dir_output_base = "/mnt/redwood/local_drive/chanud/whatchem_2.1_CNY_4cast/output_data";                             # Directory for output data

#### END USER SPECIFICATIONS ################################################################




### RUN MODEL ###############################################################################

### TEST CODE ###############
# This code is for testing a single WHATCH'EM run
if (0) {
	# Input data parameters
	my $dataSrc    = "GEOS";             # Climate data source ("MERRA_IMERG" or "GEOS")
	my $lead_time  = 0;                  # only relevant for GEOS data; gives lead time in months
	my $ens_num    = 1;                  # only relevant for GEOS data; gives ensemble number
	my $loc        = "NuwaraEliya";      # Geographic location/city name
	my $tspan_strt = "2020120100";       # Start time for model run (YYYYMMDDHH)
	my $tspan_end  = "2020123123";       # End time for model run (YYYYMMDDHH)
	my $tspan      = "${tspan_strt}_${tspan_end}";
	my $dir_input  = File::Spec->catdir($dir_input_base, $dataSrc);
	my $dir_output = File::Spec->catdir($dir_output_base, $dataSrc);
	whatchem_main($dataSrc, $loc, $tspan, $name_cont, $name_refl, $name_shade, 
				$wh_mul, $doClim, $cflag, $sflag, $mfill, $fill_intv, $mempty, $empty_intv,
				$tStep_model, $tStep_input, $missVal, $path_exe, $dir_input, $dir_output,
				$lead_time, $ens_num);
	exit;
}
### TEST CODE ###############



# Input data parameters
my @dataSrcs   = qw/ GEOS /; #qw/ MERRA_IMERG GEOS /;
my @lead_times = (0);
my @ens_nums   = (1..5);
my @locs       = qw/ Negombo Jaffna NuwaraEliya /;
my @years      = (2001..2020);
my @months     = (1..12);

# Run model for each set of input data parameters
for my $dataSrc (@dataSrcs){
	my $dir_input  = File::Spec->catdir($dir_input_base, $dataSrc);
	my $dir_output = File::Spec->catdir($dir_output_base, $dataSrc);
	for my $loc (@locs){
		for my $year (@years) {
			for my $month (@months) {
				print "Data: $dataSrc  Location: ${loc}  Year: ${year}  Month: ${month}\n";

				my $tspan = make_timespan($year, $month); # generate timespan string based on year and month of interest
				print "Timespan: ${tspan}\n";

				# Run the WHATCH'EM code
				# (For GEOS data we do this for each lead time and ensemble number)
				if ($dataSrc eq 'GEOS'){
					for my $lead_time (@lead_times){
						print " Lead time: ${lead_time}\n";
						for my $ens_num (@ens_nums){
							print " ens${ens_num}\n";
							whatchem_main($dataSrc, $loc, $tspan, $name_cont, $name_refl, $name_shade, 
								          $wh_mul, $doClim, $cflag, $sflag, $mfill, $fill_intv, $mempty, $empty_intv,
								          $tStep_model, $tStep_input, $missVal, $path_exe, $dir_input, $dir_output,
								          $lead_time, $ens_num);
						}
					}
				} else {
					whatchem_main($dataSrc, $loc, $tspan, $name_cont, $name_refl, $name_shade, 
								  $wh_mul, $doClim, $cflag, $sflag, $mfill, $fill_intv, $mempty, $empty_intv,
								  $tStep_model, $tStep_input, $missVal, $path_exe, $dir_input, $dir_output);
				}
				
			
			}

		}
	}
}


### Helper functions for model-running loop #####################

### Check whether a given year is a leap year
#
# Arguments
#   year - year to test for leap-year-ness
sub is_leap_year {
	my ($year) = @_;
	return 0 if $year % 4;         # Not divisible by 4 -> not a leap year
	return 1 if $year % 400 == 0;  # Divisible by 400 -> leap year
	return 0 if $year % 100 == 0;  # Divisible by 100 (but not by 400) -> not a leap year
	return 1;                      # Else -> leap year
}

### Check number of days in a given year + month
#
# Arguments
#   year  - year in question
#   month - month in question
sub get_days_in_month {
	my ($year, $month) = @_;
	return 31 if $month =~ /^(1|3|5|7|8|10|12)$/;
	return 30 if $month =~ /^(4|6|9|11)$/;
	return is_leap_year($year) ? 29 : 28; # Feb has 29 days in leap year, 28 days otherwise
}

### Create timespan string indicating the first to last days 
#   of a given year + month (format: YYYYMMDD_YYYYMMDD)
#
# Arguments
#   year  - year in question
#   month - month in question
sub make_timespan {
	my ($year, $month) = @_;
	my $days = get_days_in_month($year, $month);
	my $start = sprintf("%04d%02d0100", $year, $month);
	my $end   = sprintf("%04d%02d%02d23", $year, $month, $days);
	return "${start}_${end}";
}

### END RUN MODEL ##############################################################################




### Main function for running WHATCH'EM.
#
# Arguments
#       dataSrc - Source of climate data
#           loc - Geographic location/city name
#         tspan - Timespan of input data and of model run (YYYYMMDDHH_YYYYMMDDHH)
#     name_cont - Name of container
#     name_refl - Name of reflectivity (i.e., color)
#    name_shade - Name of shade condition
#        wh_mul - Multiplier of initial (water height / container height) ratio
#        doClim - Whether input data corresponds to climatology (T/F)
#         cflag - Whether cloud fraction data is input (T/F)
#         sflag - Whether soil temperature data is input (T/F)
#         mfill - Manual fill flag (T/F)
#     fill_intv - Manual fill interval (hours)
#        mempty - Manual empty flag (T/F)
#    empty_intv - Manual empty interval (hours)
#   tStep_model - Model time step (seconds)
#   tStep_input - Time interval of input data (seconds)
#       missVal - Missing value of input data
#      path_exe - Path to container.exe file
#     dir_input - Directory of input data
#  dir_output - Directory for output data
#     lead_time - indicates GEOS lead time in months (optional, only use if dataSrc=GEOS)
#       ens_num - indicates GEOS ensemble number (optional, only use if dataSrc=GEOS)
sub whatchem_main {
	my ($dataSrc, $loc, $tspan, $name_cont, $name_refl, $name_shade, 
	    $wh_mul, $doClim, $cflag, $sflag, $mfill, $fill_intv, $mempty, $empty_intv,
		$tStep_model, $tStep_input, $missVal, $path_exe, $dir_input, $dir_output,
		$lead_time, $ens_num) = @_;
	
	# If not using the GEOS-only arguments, define them as undefined.
	$lead_time = undef unless defined $lead_time;
	$ens_num   = undef unless defined $ens_num;

	# Clouds (only used if not provided)
	my @clouds = ([0.0,0.0,0.0,0.0,0.0,0.0], [0.33,0.33,0.33,0.1,0.1,0.1], [0.66,0.66,0.66,0.33,0.33,0.33]); # Inner index is low, middle, and high clouds for day, then same for night.
	my @names_cloud = qw/ Clear PartlyCloudy Overcast /; # Names of cloud cover conditions (if specified)


	#### PROCESSING OF USER SPECIFICATIONS ####
	# Create assorted variables based on the user specifications.

	# Run names
	# name indicates data source, location, and timespan of input data (plus lead time and ensemble number, if relevant)
	my $name_srcLocTime;
	if ($dataSrc eq 'GEOS') {
        $name_srcLocTime = "${dataSrc}_lead${lead_time}_ens${ens_num}-${loc}-${tspan}"; 
	} else {
		$name_srcLocTime = "${dataSrc}-${loc}-${tspan}";
	}
	
	if ($doClim eq 'T') {           # if input data is climatology
		$name_srcLocTime = "$name_srcLocTime-clim";
	}

	my $name_intv = "intv";                       # name indicating interventions
	if ($mfill eq 'T') {
		$name_intv = "${name_intv}F";
	} else {
		$name_intv = "${name_intv}0";
	}
	if ($mempty eq 'T') {
		$name_intv = "${name_intv}E";
	} else {
		$name_intv = "${name_intv}0";
	}
	
	my $name_contSpecs = "$name_cont-$name_refl-$name_shade"; # name indicating container specifications

	# Filepaths
	my $input_file = "$dir_input/container_input_$name_srcLocTime.txt";
	my $path_out_data = "$dir_output/$name_srcLocTime";
	my $output_file_base = "container_output_$name_srcLocTime";
	if (!(-d "$path_out_data")) { system "mkdir $path_out_data" }
	
	## Look up variable values based on the user inputs (see the lookup_vars() subroutine).
	#  Note: these float variables are all converted to strings via sprintf() so that they
	#        have consistent formatting. Otherwise the Fortran code may misinterpret them.
	my %all_vars = lookup_vars($loc, $name_cont, $name_refl, $name_shade);
	my %loc_vars        = %{$all_vars{"loc_vars"}};
	my %cont_shape_vars = %{$all_vars{"cont_shape_vars"}};
	my %cont_refl_vars  = %{$all_vars{"cont_refl_vars"}};
	my %cont_shade_vars = %{$all_vars{"cont_shade_vars"}};

	# Geographic location
	my $lat       = sprintf("%.2f", $loc_vars{"lat"});              # Site latitude (deg N)
	my $lon       = sprintf("%.2f", $loc_vars{"lon"});              # Site longitude (deg E)
	my $elev      = sprintf("%.2f", $loc_vars{"elev"});             # Site elevation (m)

	# Containers
	my $cshape    = $cont_shape_vars{"cshape"};                     # Shape of container
	my $radius1_t = sprintf("%.4f", $cont_shape_vars{"radius1_t"}); # Top container length 1 (in m)
	my $radius2_t = sprintf("%.4f", $cont_shape_vars{"radius2_t"}); # Top container length 2 (in m)
	my $radius1_b = sprintf("%.4f", $cont_shape_vars{"radius1_b"}); # Body container length 1 (in m)
	my $radius2_b = sprintf("%.4f", $cont_shape_vars{"radius2_b"}); # Body container length 2 (in m)
	my $height    = sprintf("%.4f", $cont_shape_vars{"height"});    # Container height (in m)
	my $conduc    = sprintf("%.2f", $cont_shape_vars{"conduc"});    # Thermal conductivity of container (in W m-1 K-1)
	my $thick     = sprintf("%.4f", $cont_shape_vars{"thick"});     # Thickness of container (in m)

	# Reflectivities (colors)
	my $refl      = sprintf("%.2f", $cont_refl_vars{"refl"});       # Reflectivity coefficient of container

	# Shade
	my $shade     = sprintf("%.2f", $cont_shade_vars{"shade"});     # Fraction of shade

	# Misc.
	my $wh      = sprintf("%.4f", $height*$wh_mul);                 # Container water height (same units as container height)
	my $t_diff  = sprintf("%.2f", $tStep_input/$tStep_model);       # Number of input data timesteps in one model timestep
	$tStep_model = sprintf("%.2f", $tStep_model);                   # Now that t_diff is calculated, convert to formatted string

	# Water and ground conditions
	# (estimate as the initial air temperature from the climate input data)
	my $init_airTemp = sprintf("%.2f", get_init_airTemp($input_file)); # get first air temperature value from input climate data (deg C)
    my $tw = $init_airTemp;         # Initial water temperature (deg C)
    my $tg = $init_airTemp;         # Initial ground temperature (deg C) (only used if sflag is FALSE, i.e. if soil temps not provided)

	#### END PROCESSING OF USER SPECIFICATIONS ####


	# Welcome message
	print "##############################\n";
	print "#### WELCOME TO WHATCH'EM ####\n";
	print "##############################\n";
	print "# Starting...\n";
	print "# Run name:       $name_srcLocTime\n";
	print "# Interventions:  $name_intv\n";
	print "# Input file:     $input_file\n";

	# Get array sizes for use throughout program
	my $exp_cloud_size = scalar (@clouds);

	# Open master statistics files
	open(MTOUTFILE, ">stats/stats_temperature.txt") or die "Cannot open master temperature stats file: $!";
	open(MWOUTFILE, ">stats/stats_waterheight.txt") or die "Cannot open master water height stats file: $!";
	if ($cflag eq 'F') {
		printf MTOUTFILE "%-30s %-10s %-15s %-10s %-15s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-6s\n", "LOC-TSPAN", "CONTAINER", "REFL", "SHADE", "CLOUDS", "TMEAN", "TMAX", "TMIN", "CMEAN", "CMAX", "CMIN", "TAMEAN", "EVAP";
		printf MWOUTFILE "%-30s %-10s %-15s %-10s %-15s %-15s\n", "LOC-TSPAN", "CONTAINER", "REFL", "SHADE", "CLOUDS", "WATER HEIGHT";
	}
	else {
		printf MTOUTFILE "%-30s %-10s %-15s %-10s %-6s %-6s %-6s %-6s %-6s %-6s %-6s %-6s\n", "LOC-TSPAN", "CONTAINER", "REFL", "SHADE", "TMEAN", "TMAX", "TMIN", "CMEAN", "CMAX", "CMIN", "TAMEAN", "EVAP";
		printf MWOUTFILE "%-30s %-10s %-15s %-10s %-15s\n", "LOC-TSPAN", "CONTAINER", "REFL", "SHADE", "WATER HEIGHT";
	}

	# Loop through each cloud (if not provided)
	if ($cflag eq 'F') {
		for (my $l = 0; $l < $exp_cloud_size; $l++) {
			my $name_cloud = $names_cloud[$l];
			print "##### $name_cloud\n";
			my $output_file = $output_file_base."-$name_intv-$name_contSpecs-$name_cloud.txt";
			my @clouds_exp = @{$clouds[$l]};
			my $cld = sprintf("%.2f", $clouds_exp[0]);
			my $cmd = sprintf("%.2f", $clouds_exp[1]);
			my $chd = sprintf("%.2f", $clouds_exp[2]);
			my $cln = sprintf("%.2f", $clouds_exp[3]);
			my $cmn = sprintf("%.2f", $clouds_exp[4]);
			my $chn = sprintf("%.2f", $clouds_exp[5]);

			# Construct parameter file
			print "##### Constructing parameter file...\n";
			open OUTFILE, ">param_file.txt";
			print OUTFILE "$input_file\n";
			print OUTFILE "$loc\n";
			print OUTFILE "$tspan\n";
			print OUTFILE "$lat\n";
			print OUTFILE "$lon\n";
			print OUTFILE "$elev\n";
			print OUTFILE "$cshape\n";
			print OUTFILE "$radius1_t\n";
			print OUTFILE "$radius2_t\n";
			print OUTFILE "$radius1_b\n";
			print OUTFILE "$radius2_b\n";
			print OUTFILE "$height\n";
			print OUTFILE "$refl\n";
			print OUTFILE "$conduc\n";
			print OUTFILE "$thick\n";
			print OUTFILE "$wh\n";
			print OUTFILE "$tw\n";
			print OUTFILE "$tg\n";
			print OUTFILE "$tStep_model\n";
			print OUTFILE "$t_diff\n";
			print OUTFILE "$shade\n";
			print OUTFILE "$cflag\n";
			print OUTFILE "$sflag\n";
			print OUTFILE "$mfill\n";
			print OUTFILE "$fill_intv\n";
			print OUTFILE "$mempty\n";       # (ADDED ON 2023-08-09, CNY)
			print OUTFILE "$empty_intv\n";   # (ADDED ON 2023-08-09, CNY)
			print OUTFILE "$output_file\n";
			print OUTFILE "$cld\n";
			print OUTFILE "$cmd\n";
			print OUTFILE "$chd\n";
			print OUTFILE "$cln\n";
			print OUTFILE "$cmn\n";
			print OUTFILE "$chn\n";
			close OUTFILE;

			# Run container modeling program
			print "##### Running container model...";
			system "$path_exe/container.exe";

			# Move output file to correct directory
			print "##### Writing output files...\n";
			system "mv $output_file $path_out_data";

			# Run daily energy balance and basic stats
			my $einfile = "$path_out_data/$output_file";
			my $eoutfile = "$path_out_data/energybal_avg_$name_srcLocTime-$name_intv-$name_contSpecs-$name_cloud.txt";
			my $soutfile = "$path_out_data/stats_$name_srcLocTime-$name_intv-$name_contSpecs-$name_cloud.txt";
			&daily_energy_balance_stats($einfile, $eoutfile, $soutfile, $name_srcLocTime, $name_cont, $name_refl, $name_shade, $cflag, $name_cloud);
		}
	} else {
		my $output_file = $output_file_base."-$name_intv-$name_contSpecs.txt";
		# Construct parameter file
		print "#### Constructing parameter file...\n";
		open OUTFILE, ">param_file.txt";
		print OUTFILE "$input_file\n";
		print OUTFILE "$loc\n";
		print OUTFILE "$tspan\n";
		print OUTFILE "$lat\n";
		print OUTFILE "$lon\n";
		print OUTFILE "$elev\n";
		print OUTFILE "$cshape\n";
		print OUTFILE "$radius1_t\n";
		print OUTFILE "$radius2_t\n";
		print OUTFILE "$radius1_b\n";
		print OUTFILE "$radius2_b\n";
		print OUTFILE "$height\n";
		print OUTFILE "$refl\n";
		print OUTFILE "$conduc\n";
		print OUTFILE "$thick\n";
		print OUTFILE "$wh\n";
		print OUTFILE "$tw\n";
		print OUTFILE "$tg\n";
		print OUTFILE "$tStep_model\n";
		print OUTFILE "$t_diff\n";
		print OUTFILE "$shade\n";
		print OUTFILE "$cflag\n";
		print OUTFILE "$sflag\n";
		print OUTFILE "$mfill\n";
		print OUTFILE "$fill_intv\n";
		print OUTFILE "$mempty\n";       # (ADDED ON 2023-08-09, CNY)
		print OUTFILE "$empty_intv\n";   # (ADDED ON 2023-08-09, CNY)
		print OUTFILE "$output_file\n";
		close OUTFILE;

		# Run container modeling program
		print "#### Running container model...";
		system "$path_exe/container.exe";

		# Move output file to correct directory
		print "#### Writing output files...\n";
		system "mv $output_file $path_out_data";

		# Run daily energy balance and basic stats
		my $einfile = "$path_out_data/$output_file";
		my $eoutfile = "$path_out_data/energybal_avg_$name_srcLocTime-$name_intv-$name_contSpecs.txt";
		my $soutfile = "$path_out_data/stats_$name_srcLocTime-$name_intv-$name_contSpecs.txt";
		&daily_energy_balance_stats($einfile, $eoutfile, $soutfile, $name_srcLocTime, $name_cont, $name_refl, $name_shade, $cflag);
	}

	close MTOUTFILE;
	close MWOUTFILE;

	# Now open the summary statistics files to run some additional plots
	print "# Writing summary statistics files and plots...\n";
	my %data;
	my @tvals = qw/ TMEAN TMAX TMIN CMEAN CMAX CMIN /;
	open MTINFILE, "<stats/stats_temperature.txt" or die "Cannot open input stats temperature file: $!";
	open MWINFILE, "<stats/stats_waterheight.txt" or die "Cannot open input stats water height file: $!";
	my $junk = <MTINFILE>;
	$junk = <MWINFILE>;
	while (my $tline = <MTINFILE>) {
		my $wline = <MWINFILE>;
		my @tparts = split /\s+/, $tline; # split at whitespace
		my @wparts = split /\s+/, $wline; # split at whitespace
		if ($cflag eq 'F') {
			$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{$tparts[4]}->{'TMEAN'} = $tparts[5];
			$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{$tparts[4]}->{'TMAX'} = $tparts[6];
			$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{$tparts[4]}->{'TMIN'} = $tparts[7];
			$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{$tparts[4]}->{'CMEAN'} = $tparts[8];
			$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{$tparts[4]}->{'CMAX'} = $tparts[9];
			$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{$tparts[4]}->{'CMIN'} = $tparts[10];
			$data{$wparts[0]}->{$wparts[1]}->{$wparts[2]}->{$wparts[3]}->{$tparts[4]}->{'WH'} = $wparts[5];
		} else {
			$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{'TMEAN'} = $tparts[4];
			$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{'TMAX'} = $tparts[5];
			$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{'TMIN'} = $tparts[6];
			$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{'CMEAN'} = $tparts[7];
			$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{'CMAX'} = $tparts[8];
			$data{$tparts[0]}->{$tparts[1]}->{$tparts[2]}->{$tparts[3]}->{'CMIN'} = $tparts[9];
			$data{$wparts[0]}->{$wparts[1]}->{$wparts[2]}->{$wparts[3]}->{'WH'} = $wparts[4];
		}
	}
	close MTINFILE;
	close MWINFILE;

	open WTEMPFILE, ">temp_plot_wh.txt" or die "Cannot open temp plot file: $!";
	foreach my $tval (@tvals) {
		open TTEMPFILE, ">temp_plot_t_$tval.txt" or die "Cannot open temp plot file: $!";
		if ($cflag eq 'F') {
			foreach my $name_cloud (@names_cloud) {
				printf TTEMPFILE "%6.2f ", $data{$name_srcLocTime}->{$name_cont}->{$name_refl}->{$name_shade}->{$name_cloud}->{$tval};
			}
		} else {
			printf TTEMPFILE "%6.2f ", $data{$name_srcLocTime}->{$name_cont}->{$name_refl}->{$name_shade}->{$tval};
		}
		print TTEMPFILE "\n";
		close TTEMPFILE;
	}

	if ($cflag eq 'F') {
		foreach my $name_cloud (@names_cloud) {
			printf WTEMPFILE "%6.2f ", $data{$name_srcLocTime}->{$name_cont}->{$name_refl}->{$name_shade}->{$name_cloud}->{'WH'};
		}
	} else {
		printf WTEMPFILE "%6.2f ", $data{$name_srcLocTime}->{$name_cont}->{$name_refl}->{$name_shade}->{'WH'};
	}

	print WTEMPFILE "\n";
	close WTEMPFILE;

	# Clean up
	print "# Cleaning up...\n";
	system "rm temp_plot_t_*.txt temp_plot_wh.txt param_file.txt";
	print "# Finished!\n";

}







### Helper functions for whatchem_main() #####################

### lookup_vars()
#     Look up specific data values based on the more general user inputs
#     (e.g., get lat and lon given the city name).
#  
sub lookup_vars {
	my ($loc, $name_cont, $name_refl, $name_shade) = @_;

	# Create hashes containing data corresponding to all possible user inputs.
	my %loc_lookupData = (
		"Negombo"     => {"lat" => 7.2008, "lon" => 79.8737, "elev" => 2.},
		"Jaffna"      => {"lat" => 9.6615, "lon" => 80.0255, "elev" => 5.},
		"NuwaraEliya" => {"lat" => 6.9497, "lon" => 80.7891, "elev" => 1868.},
		"Trincomalee" => {"lat" => 8.5874, "lon" => 81.2152, "elev" => 8.}
	);
	my %cont_shape_lookupData = (
		"Coffee" => {"cshape" => "ROUND", "radius1_t" => 0.082, "radius2_t" => 0.082,
					"radius1_b" => 0.082, "radius2_b" => 0.082, "height" => 0.191,
					"conduc" => 0.30, "thick" => 0.0018},
		"Bucket" => {"cshape" => "ROUND", "radius1_t" => 0.131, "radius2_t" => 0.131,
					"radius1_b" => 0.131, "radius2_b" => 0.131, "height" => 0.368,
					"conduc" => 0.50, "thick" => 0.0023},
		"Barrel" => {"cshape" => "ROUND", "radius1_t" => 0.298, "radius2_t" => 0.298,
					"radius1_b" => 0.298, "radius2_b" => 0.298, "height" => 0.921,
					"conduc" => 0.50, "thick" => 0.0048}
	);
	my %cont_refl_lookupData = (
		"Black" => {"refl" => 0.1},
		"Gray"  => {"refl" => 0.5},
		"White" => {"refl" => 0.9}
	);
	my %cont_shade_lookupData = (
		"NoShade"   => {"shade" => 0.00},
		"HalfShade" => {"shade" => 0.50},
		"FullShade" => {"shade" => 1.00}
	);

	# Look up and return the specific values corresponding to the actual user inputs.
	my %loc_vars        = %{$loc_lookupData{$loc}};
	my %cont_shape_vars = %{$cont_shape_lookupData{$name_cont}};
	my %cont_refl_vars  = %{$cont_refl_lookupData{$name_refl}};
	my %cont_shade_vars = %{$cont_shade_lookupData{$name_shade}};
    
	my %all_vars = ( # store /references/ to hashes here (trying to store hashes directly doesn't work)
		"loc_vars"        => \%loc_vars,
		"cont_shape_vars" => \%cont_shape_vars,
		"cont_refl_vars"  => \%cont_refl_vars,
		"cont_shade_vars" => \%cont_shade_vars
	);
	return %all_vars;
}




### Get initial air temperature from input climate data file.
#   Note: Assumes the name of the air temperature column, which is
#         hardcoded within the $varName_airTemp variable.
#
# Arguments
#   input_file - input file of climate data
sub get_init_airTemp {
	my ($input_file) = @_;
	
	# read relevant lines from input file
    open(FH, '<', $input_file) or die "Error on file open: $!"; # open input file
    my $line_varNames         = <FH>; # read first line, which contains the variable names as column headings
    my $line_varVals_firstRow = <FH>; # read second line, which is first line of data
    close(FH);

    # create hash of the variable names and their values
	# (this simplifies finding the correct column in the next step)
    my @varNames         = split(" ", $line_varNames);
    my @varVals_firstRow = split(" ", $line_varVals_firstRow);
    my %inputData_firstRow;
    @inputData_firstRow{@varNames} = @varVals_firstRow;

    # get initial air temperature value
	my $varName_airTemp = "T"; # this must match the name of the air temperature column
    my $init_airTemp = $inputData_firstRow{$varName_airTemp};

	return $init_airTemp;
}




### daily_energy_balance_stats()
#     Calculate energy balance statistics.
#    
sub daily_energy_balance_stats {
	my ($einfile, $eoutfile, $soutfile, $name_srcLocTime, $name_cont, $name_refl, $name_shade, $cflag, $name_cloud) = @_;
	# Construct average daily energy balance and some other statistics
	my %ebal;
	my %tw;
	my %tc;
	my %ta;
	my %evap;
	my $whc = 0;
	open EINFILE, "<$einfile" or die "Cannot open input file for energy balance: $!";
	for (my $x = 0; $x < 7; $x++) {
		my $junk = <EINFILE>; # Header lines
	}
	while (my $line = <EINFILE>) {
		chomp ($line);
		my @parts = split /\s+/, $line;
		my $md = $parts[1].$parts[2];
		my $hour = $parts[3];
		if ($parts[24] != $missVal) {
			push (@{$ebal{$hour}->{'SWDOWN'}}, $parts[4]);
			push (@{$ebal{$hour}->{'SWSIDE'}}, $parts[5]);
			push (@{$ebal{$hour}->{'LWDOWN'}}, $parts[6]);
			push (@{$ebal{$hour}->{'LWUP'}}, $parts[7]);
			push (@{$ebal{$hour}->{'LWIN'}}, $parts[8]);
			push (@{$ebal{$hour}->{'LWOUT'}}, $parts[9]);
			push (@{$ebal{$hour}->{'SHUP'}}, $parts[10]);
			push (@{$ebal{$hour}->{'SHOUT'}}, $parts[11]);
			push (@{$ebal{$hour}->{'LHUP'}}, $parts[12]);
			push (@{$ebal{$hour}->{'GHW'}}, $parts[13]);
			push (@{$ebal{$hour}->{'GHC'}}, $parts[14]);
			push (@{$ebal{$hour}->{'COND'}}, $parts[15]);
			push (@{$ebal{$hour}->{'BAL_W'}}, $parts[16]);
			push (@{$ebal{$hour}->{'BAL_C'}}, $parts[17]);
			my $sw = $parts[4] + $parts[5];
			push (@{$ebal{$hour}->{'SW'}}, $sw);
			my $lw = $parts[6] - $parts[7] + $parts[8] - $parts[9];
			push (@{$ebal{$hour}->{'LW'}}, $lw);
			my $sh = -1.*($parts[10] + $parts[11]);
			push (@{$ebal{$hour}->{'SH'}}, $sh);
			push (@{$evap{$md}}, $parts[19]);
			push (@{$ebal{$hour}->{'TG'}}, $parts[21]);
			push (@{$ebal{$hour}->{'TA'}}, $parts[22]);
			push (@{$ebal{$hour}->{'TCON'}}, $parts[23]);
			push (@{$ebal{$hour}->{'TW'}}, $parts[24]);
			push (@{$tw{$md}}, $parts[24]);
			push (@{$tc{$md}}, $parts[23]);
			push (@{$ta{$md}}, $parts[22]);
			if ($parts[20] >= 0.015) {
				$whc++;
			}
		}
	}
	close EINFILE;

	# Write out daily average energy balance file
	open EOUTFILE, ">$eoutfile" or die "Cannot open output file for energy balance: $!";
	print EOUTFILE "Hour,SWDOWN,SWSIDE,LWDOWN,LWUP,LWIN,LWOUT,SHUP,SHOUT,LHUP,GHW,GHC,COND,BAL_W,BAL_C,SW,LW,SH,TG,TA,TCON,TW\n";
	for (my $hr = 0; $hr < 24; $hr++) {
		printf EOUTFILE "%02s ", $hr;
		my @terms = qw/ SWDOWN SWSIDE LWDOWN LWUP LWIN LWOUT SHUP SHOUT LHUP GHW GHC COND BAL_W BAL_C SW LW SH TG TA TCON TW /;
		foreach my $term (@terms) {
			my $estats = Statistics::Descriptive::Full->new();
			$estats->add_data(@{$ebal{$hr}->{$term}});
			my $mean = $estats->mean();
			printf EOUTFILE "%8.2f ", $mean;
		}
		print EOUTFILE "\n";
	}
	close EOUTFILE;

	# Calculate some basic statistics and output to file
	open SOUTFILE, ">$soutfile" or die "Cannot open output file for stats: $!";
	my %tastats;
	$tastats{'TMEAN'} = Statistics::Descriptive::Full->new();
	$tastats{'TMAX'} = Statistics::Descriptive::Full->new();
	$tastats{'TMIN'} = Statistics::Descriptive::Full->new();
	$tastats{'CMEAN'} = Statistics::Descriptive::Full->new();
    $tastats{'CMAX'} = Statistics::Descriptive::Full->new();
    $tastats{'CMIN'} = Statistics::Descriptive::Full->new();
	$tastats{'TAMEAN'} = Statistics::Descriptive::Full->new();
	$tastats{'EVAP'} = Statistics::Descriptive::Full->new();

	foreach my $md (sort keys %tw) {
		my $tstats = Statistics::Descriptive::Full->new();
		my $astats = Statistics::Descriptive::Full->new();
		my $cstats = Statistics::Descriptive::Full->new();
		$tstats->add_data(@{$tw{$md}});
		$astats->add_data(@{$ta{$md}});
		$cstats->add_data(@{$tc{$md}});
		$tastats{'EVAP'}->add_data(@{$evap{$md}});
		my $tmean = $tstats->mean();
		my $tmax = $tstats->max();
		my $tmin = $tstats->min();
		my $tamean = $astats->mean();
		my $cmean = $cstats->mean();
        my $cmax = $cstats->max();
        my $cmin = $cstats->min();
		if (scalar(@{$tw{$md}}) == 24) {
			$tastats{'TMEAN'}->add_data($tmean);
			$tastats{'TMAX'}->add_data($tmax);
			$tastats{'TMIN'}->add_data($tmin);
			$tastats{'CMEAN'}->add_data($cmean);
            $tastats{'CMAX'}->add_data($cmax);
            $tastats{'CMIN'}->add_data($cmin);
			$tastats{'TAMEAN'}->add_data($tamean);
		}
	}
	print SOUTFILE "Average Water Daily Temperature (Mean, Max, Min): \n";
	my $tmean = $tastats{'TMEAN'}->mean();
	my $tmax = $tastats{'TMAX'}->mean();
	my $tmin = $tastats{'TMIN'}->mean();
	my $cmean = $tastats{'CMEAN'}->mean();
    my $cmax = $tastats{'CMAX'}->mean();
    my $cmin = $tastats{'CMIN'}->mean();
	my $tamean = $tastats{'TAMEAN'}->mean();
	my $evapmean = $tastats{'EVAP'}->mean();
	printf SOUTFILE "%6.2f %6.2f %6.2f\n", $tmean, $tmax, $tmin;
	print SOUTFILE "Average Container Daily Temperature (Mean, Max, Min): \n";
	printf SOUTFILE "%6.2f %6.2f %6.2f\n", $cmean, $cmax, $cmin;
	if ($cflag eq 'F') {
		printf MTOUTFILE "%-50s %-10s %-15s %-10s %-15s %-6.2f %-6.2f %-6.2f %-6.2f%-6.2f %-6.2f %-6.2f %-6.4f\n", $name_srcLocTime, $name_cont, $name_refl, $name_shade, $name_cloud, $tmean, $tmax, $tmin, $cmean, $cmax, $cmin, $tamean, $evapmean;
	} else {
		printf MTOUTFILE "%-50s %-10s %-15s %-10s %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f %-6.2f %-6.4f\n", $name_srcLocTime, $name_cont, $name_refl, $name_shade, $tmean, $tmax, $tmin, $cmean, $cmax, $cmin, $tamean, $evapmean;
	}
	print SOUTFILE "Number of days water level > 15 mm: \n";
	$whc = $whc / 24;
	printf SOUTFILE "%6.2f\n", $whc;
	if ($cflag eq 'F') {
		printf MWOUTFILE "%-50s %-10s %-15s %-10s %-15s %-6.2f\n", $name_srcLocTime, $name_cont, $name_refl, $name_shade, $name_cloud, $whc;
	} else {
		printf MWOUTFILE "%-50s %-10s %-15s %-10s %-6.2f\n", $name_srcLocTime, $name_cont, $name_refl, $name_shade, $whc;
	}
	print SOUTFILE "Frequency Distribution: \n";
	my @bins = (10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50); ### added more bins (2023-??-??, CNY)
	my $fd = $tastats{'TMEAN'}->frequency_distribution_ref(\@bins);
	for (sort {$a <=> $b} keys %{$fd}) {
		print SOUTFILE "$_  $fd->{$_}\n";
	}
	close SOUTFILE;
}
