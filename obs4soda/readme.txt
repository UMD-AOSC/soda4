nodc2soda_55levs.f     -  this subroutine reads WODB's temperature and salinty observation profile date 
                          for the depth 0-2000 meters, interpolates it to the SODA4 vertical levels, 
		          and writes it in ASCII text files.

t_dat2qc_tpot_55levs.f - this subroutine performes quality control for the potential temperature profiles.

s_dat2qc_55levs_2.f    - this subroutine performes quality control for the salinity profiles.

qc2bin_55levs.f        - this subroutine bin profiles by 1degree x 1degree x 5days bins
