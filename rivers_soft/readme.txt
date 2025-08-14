gridded_discgarge_720x360x528.m         - this is main MATLAB subroutine to create 1900-2023 years 
                                          0.5x0.5 degree global river discharge data set
				          (gridded_discharge_720x360x528.nc)

Dai_Tr_extension_filled_adjusted_2023.m	- this is MATLAB subroutine to:
                                          1. Update Dai_Trenberth data file by 20 new river observations
					     from CRDC
					  2. Add Arctic rivers observations 
                                          3. Adjust rivers discharge to the unmonitored areas in each 
					     ocean follow to Dai and Trenberth (2002)

read_Arctic_river.m                     - this is MATLAB subroutine to read Arctic data from separate
                                          text files
					  
read_crdc_river.m                       - this is MATLAB subroutine to read CRDC data from separate 
                                          EXCEL's text files


