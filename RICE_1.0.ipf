#pragma TextEncoding = "Windows-1252"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma rtGlobals=3		// Use modern global access method and strict wave access.



//The functions described below comprise of the program described as the Refractive Index Confidence Explorer or RICE.  RICE relies on an interative approach
//to estimate uncertainty in optical closure calculations using discrete solutions.  For each RICE cycle RICE attempts to determine if a given discrete refractive index value could have provided
//the refractive index that was calculated experimentally( or test values).  By doing this a large number of times for a series of near refractive index values, 
// a probability distribution is created of refractive index values that could have yielded the observed value under the uncertainty conditions.  All calculations operate under
// assumptions of Mie Theory.

//RICE was developed at the University of California Riverside in 2020 by Alexander Frie and Roya Bahreini.  RICE was developed with support from the national science foundation (AGS 1454374).
//If you have questions feel free to contact Alex at afrie003@ucr.edu.  


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

Menu "RICE"
	"RICE 1.0", User_interface_RICE()
	"Refractive Index Calcuation",  RI_calculation() 
	"Reanalyze RICE Interval", RICE_Reanalyze()
	//creates a menu item that iniates the macro
End


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


function RICE_Reanalyze()

string adaptive_output="RICE_Output"
variable CI=0.95

prompt adaptive_output, "RICE_Output Matrix, should be named 'RICE_Output'", popup, WaveList("*",";","")
prompt CI, "New confidence level (0-1)"
doprompt "RICE_Reanalyze", adaptive_output, CI

if(v_flag==1) //aborts on cancel
	print "canceled by user"
	abort 
endif



Output_updater_k($Adaptive_output,CI)
Output_updater_n($Adaptive_output,CI)







end


function User_interface_RICE() //Creates a simple user interface for the selection of input waves

string dia_bins="dia_bins"  
string Size_t_series="Size_t_series"
string n_t_series="RI_t_series"
string k_t_series="k_t_series"



prompt dia_bins, "Particle Diameter Bin Locations", popup, WaveList("*",";","") 
prompt Size_t_series, "2D size Distribtuion (d(N)/dlog(d))--Columns should be time and Rows should be particle size.", popup, WaveList("*",";","") 
prompt n_t_series,"n value time series", popup, WaveList("*",";","") 
prompt k_t_series,"k value time series", popup, WaveList("*",";","")
	
doprompt "Input Data", dia_bins, Size_t_series, n_t_series, k_t_series // prompts for the input data

if(v_flag==1) //aborts on cancel
	print "canceled by user"
	abort 
endif

Adapt_RI_Interval($dia_bins,$Size_t_series,$n_t_series,$k_t_series)

end
	
	
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Adapt_RI_interval prompts for the uncertainty and run conditions for the RICE analysis.  This function also repeatedly adjusts the test RI range input range in an effort to capture the probability space around the calculated refractive index.
//Adapt_RI_interval starts by inputting a small number of iterations and discrete RI values into the RI interval function.  Once the output from RI interval satisfies a set of parameters, in an effort to guess the appropriate sampling space, Adapt_RI_interval 
//increases the number of iterations and discrete RI values input into RI_interval for the final analysis of a given point.

function Adapt_RI_Interval(dia_bins, Size_t_series, RI_t_series, k_t_series)
wave dia_bins, Size_t_series, RI_t_series,K_t_series

variable ad, af=0, as=0, n_in,k_in,  k_change, k_width1, n_change, n_width1,k_step1, n_step1, ty, tx, qi, V_mid_point, iter, spin_up_cycles=5, spin_up_matrix=64, duration=0, iq, good_count=0, spin_up_repeat=1, abs_limit=0.05,scat_limit=0.05,step_modifier, matrix_final=144, SD_Dia=0.03,SD_N=0.1,SD_Scat=0.05,SD_Abs=0.05,iter_start=10, wavelength=375, iter_final=100,k_start=0, k_end=0.1, n_start=1.0, n_end=2, n_step=0.01,k_step=0.001, indicators=0, repeat_count=0,n_width=0.05, k_width=0.002, adaptors=0, bin_counter, bin_counter2,start_time
 variable idle_n, idle_k,CL=95, Calc_type=1


//promts for uncertainty conditions and wavelength of measurement

prompt SD_Dia, "Uncertainty in particle diameters"
prompt SD_N, "Uncertainty in number concentrations"
prompt SD_Scat, "Uncertainty in scattering coefficients"
prompt SD_abs, "Uncertainty in absorption coefficients"
prompt wavelength, "Wavelength of optical coefficients"
doprompt "Uncertainties" SD_Dia,SD_N,SD_Scat,SD_Abs, wavelength
if(v_flag==1) //aborts on cancel
print "canceled by user"
abort 
endif

// prompts for conditions of the analysis and output confidence interval size

prompt iter_start, "Number of iterations to perform during spin up"
prompt spin_up_matrix, "i x j value during spin up"
prompt iter_final, "Number of interations to perform when calcuating final solution"
prompt matrix_final, "i x j value when calcuating the final solution"
Prompt calc_type, "RI Calcuation type, 0=difference based (Frie et al. 2019), 1=\"chi squared based\" (Dingle et al. 2019)"
prompt CL, "Width of the output confidence interval (between 0 and 100)"
	doprompt "RICE Operational and Output Settings",iter_start, spin_up_matrix, iter_final, Matrix_final, CL 

if(v_flag==1) //aborts on cancel
print "canceled by user"
abort 
endif

CL/=100



//prompts for the name of, and creates, the output folder
killstrings /z T_name
variable run_type
string T_name="output_name"
prompt Run_type, "Do you want to keep all run files(1) or only output files(0)?  RICE outputs are much smaller if only output files are kept.    "
prompt T_name, "Output folder name"
doprompt "Output", T_name, run_type
if(v_flag==1) //aborts on cancel
print "canceled by user"
abort 
endif


killwaves/z Adaptive_output
newdatafolder /o :$T_name
setdatafolder :$T_name

make /n=(dimsize(RI_t_series,0)), space_n, space_k, k_Upper_Bound, k_Lower_Bound, n_Upper_Bound, n_Lower_Bound, count_

duplicate /o RI_t_series, n_series
duplicate /o k_t_series, k_series
duplicate /o size_t_series, Size_dist_series




///Establishes a wave to save the initial conditions so they can be repeated at each point with
variable matrix_max=spin_up_matrix
killwaves /z settings1
make /n=(1,16) settings1
settings1[0][0]=SD_Dia
settings1[0][1]=SD_N
settings1[0][2]=SD_scat
settings1[0][3]=SD_Abs
settings1[0][4]=iter_start
settings1[0][5]=n_width
settings1[0][6]=k_width
settings1[0][7]=k_step
settings1[0][8]=n_step
settings1[0][11]=spin_up_cycles
settings1[0][12]=spin_up_repeat
settings1[0][13]=abs_limit
settings1[0][14]=scat_limit
settings1[0][15]=Matrix_max




print "uncertainty conditions:", " SD Diameter=",SD_Dia," SD Number concentration=",SD_N," SD Scattering=",SD_scat," SD Absorption=",SD_Abs


variable start_=0
for (iq=0;iq<numpnts(RI_t_series);iq+=1) //performs the analysis for each point within the t_series

adaptors=0

if (RI_t_series[iq]>0) // doesn't analyze nan or negative n values
	saveexperiment
 //checks to see if this is the frist point being analyzed, if it is, it saves the settings in the settings1 wave
	if( start_==0)
		start_=1

		SD_Dia=settings1[0]
		SD_N=settings1[1]
		SD_scat=settings1[2]
		SD_Abs=settings1[3]
		iter=settings1[4]
		n_width=settings1[5]
		k_width=settings1[6]
		k_step=settings1[7]
		n_step=settings1[8]
		spin_up_cycles=settings1[11]
		spin_up_repeat=settings1[12]
		abs_limit=settings1[13]
		scat_limit=settings1[14]
		Matrix_max=settings1[15]
		
	endif

//creates a folder name for the given point and places all outputs in that folder
	killstrings /z folder_name
	string folder_name="point_"+num2str(iq)

	newdatafolder /o $folder_name
	killwaves/z settings_tracker2
	setdatafolder $folder_name


	
	
	make /o /n=1 final_run, matrix_max1
	duplicate /o settings1, settings


//sets a base value for the n and k width to explore for probable true RI values
	n_width=((SD_Dia+SD_N+SD_Scat+SD_abs)*2)*(RI_t_series[iq]-1)
	k_width=((SD_Dia+SD_N+SD_Scat+SD_abs)*2)*k_t_series[iq]


//explores the edges and center of the base n and k widths to examine which diameters are important:  removes diameters with <1% of the total absorption and scattering in any of the cases
	duplicate /o dia_bins, size_profile
	killwaves/z matrix_max1

	for(qi=0;qi<numpnts(dia_bins);qi+=1)
		size_profile[qi]=size_t_series[qi][iq]
	endfor

	bin_trim(Ri_t_series[iq],K_t_series[iq],n_width,k_width, dia_bins,size_profile,wavelength)
	wave size_input, bin_input


// begins the iterative loop where RI values near the input retrieved RI value are tested for the ability to produce the observed value under the uncertainty conditions.
//  First, the code "spins up" and establishes a likely complete sampling space using relatively quick settings then the code runs a final run with a goal of giving a confident and reproducible answer

	for (adaptors=0;adaptors<1000;adaptors+=1)


		iter=iter_start
		start_time=datetime





		if (adaptors==0)

			matrix_max=spin_up_matrix

			variable RI_points=matrix_max
		else
		
		
			
			if (RI_points<0.5*matrix_max)
			print"less 50% of n and k bins are physically meaningful, resetting to original sampling space"
				
				n_width=((SD_Dia+SD_N+SD_Scat+SD_abs)*2)*(RI_t_series[iq]-1)
				k_width=((SD_Dia+SD_N+SD_Scat+SD_abs)*2)*k_t_series[iq]

			endif
			

			
			
			
			
			variable k_mod
			if(k_t_series[iq]/0.5<1)

				k_mod=k_t_series[iq]/0.5

			else
				k_mod=1
			endif
//After one run, it checks to see if the fit parameters are satisfied, if they are, the final RI sampling space has been selected. If not, the RI sampling area is adjusted
			wave gauss_out, mid_point


			if (gauss_out[mid_point[0]][1]==-999)
			
			print"unable to fit n result with signmoid, resetting to original sampling space"
				
				n_width=((SD_Dia+SD_N+SD_Scat+SD_abs)*2)*(RI_t_series[iq]-1)
				k_width=((SD_Dia+SD_N+SD_Scat+SD_abs)*2)*k_t_series[iq]
				
			elseif (gauss_out[mid_point[0]][10]==-999)	
			print"unable to fit k result with signmoid, resetting to original sampling space"
				
				n_width=((SD_Dia+SD_N+SD_Scat+SD_abs)*2)*(RI_t_series[iq]-1)
				k_width=((SD_Dia+SD_N+SD_Scat+SD_abs)*2)*k_t_series[iq]


			endif	
			






				if (gauss_out[mid_point[0]][26]<=1)
					if (gauss_out[mid_point[0]][28]<=1)
						if (idle_n<=0.65)
							if (idle_k<=0.65)


							good_count+=1


							else 
							good_count=0
							endif
						else 
						good_count=0
						endif
					else 
					good_count=0
					endif
				else 
				good_count=0
				endif		
		

		
		
			if (good_count==1)
				
				duplicate /o dia_bins, size_profile
				killwaves/z matrix_max1

				for(qi=0;qi<numpnts(dia_bins);qi+=1)
					size_profile[qi]=size_t_series[qi][iq]
				endfor
				
				
				iter=iter_final
				adaptors=120000
				good_count*=0
				matrix_max=matrix_final
				final_run=1

				bin_trim(Ri_t_series[iq],K_t_series[iq],n_width,k_width, dia_bins,size_profile,wavelength)
				wave size_input, bin_input

			else

				final_run=0


//expands n width based on the fit of the sigmoid of the previous run results.  The max exansion is capped at 20% of the orginial value
				n_width1=n_width

	


				if (gauss_out[mid_point[0]][26]>1)
					n_width = n_width+(n_width*(gauss_out[mid_point[0]][3]-gauss_out[mid_point[0]][1]-1)*0.1)
					n_change=n_width-n_width1

					If(0.3<n_change/n_width1)
						n_width=1.2*n_width1
					endif


//adds an additional n_step of width if there was data in the last n bin of the previous run
					if (gauss_out[mid_point[0]][33]==1)
						n_width+=n_step
					endif


//if no width expansion was needed, the number of idle n bins are checked and the width is reduced if too many bins are unused


				elseif (idle_n>(0.65))

	
					n_width*=0.8

					gauss_out[mid_point[0]][29]=3


				elseif (idle_n>(0.75))
	
	

					n_width*=0.9




				endif




//changes the width of k exploration if many k values seem to lie within a single probability bin
				k_width1 = k_width
				if (gauss_out[mid_point[0]][28]>1)

					k_width = k_width+(k_width*(gauss_out[mid_point[0]][12]-gauss_out[mid_point[0]][10]-1)*0.1)
					k_change=k_width-k_width1


					If(0.3<k_change/k_width1)
						k_width=1.2*k_width1
					endif





					if (k_t_series[iq]-k_width>0)
					elseif (gauss_out[mid_point[0]][34]==1)
						k_width+=k_step
					endif



//if no width expansion was needed, the number of idle n bins are checked and the width is reduced if too many bins are unused



				elseif (idle_k>(0.65))

	
	

					k_width*=0.8
					gauss_out[mid_point[0]][30]=3



				elseif (idle_k>(0.75))
	

					k_width*=0.9
					gauss_out[mid_point[0]][30]=3
	




				endif
			endif	
		endif


//recalculates n_step and k_step based on the updates to n or k width

		n_step=(n_width/(trunc(sqrt(matrix_max)/2)))
		k_step=(k_width/(trunc(sqrt(matrix_max)/2)))



//calcuates the parameters for creating a range of test RI values.  Then creates the aformentioned range of test RI values.
	
		k_end=k_t_series[iq]+k_width
		k_start=k_t_series[iq]-k_width
		n_end=RI_t_series[iq]+n_width
		n_start=RI_t_series[iq]-n_width

		make /o /n=1 matrix_max1
		matrix_max1[0]=matrix_max

		makeRIrange(k_start,k_end,k_step,n_start,n_end,n_step)


		wave n_range, k_range

		duplicate /o n_range, n_input, k_input



// Trims the input n and k values to those that are above 1 and 0, respectively	
		bin_counter=0
		bin_counter2=-1
		for(qi=0;qi<numpnts(n_range);qi+=1)
	
			if (k_range[qi]>=0)
				if (n_range[qi]>=1)

					bin_counter=bin_counter2+1
					bin_counter2=bin_counter
					n_input[bin_counter]=n_range[qi]
					k_input[bin_counter]=k_range[qi]
				endif
			endif
  		endfor
		redimension /N=(bin_counter+1), n_input
		redimension /N=(bin_counter+1), k_input
	
	 RI_points=numpnts(n_input)
		if (numpnts(k_input)<RI_points)
			RI_points=numpnts(k_input)
		endif
	
	
	// determines the location of the midpoint value in the RI range, (the location of the input retrieved RI)
		killwaves /z mid_point
		make /n=1, mid_point
		n_in=RI_t_series[iq]
		k_in=K_t_series[iq]
		for(af=0;af<RI_points;af+=1)
			if (abs(n_in-n_input[af])<n_step/10)
				if ((abs(k_in-k_input[af])<k_step/10))
					mid_point[0]=af
					V_mid_point=af
				endif
			endif
		endfor


// print "number of iterations=",iter,"n range=",RI_t_series[iq],"±", n_width,"k range=", k_t_series[iq],"±", k_width," n values every", n_step,"k values every", k_step
		print "point",iq+1,"of", numpnts(RI_t_series)


//creates a settings wave which saves the uncertainty inputs and iteration levels for a run.
		if (adaptors==0)
			make/o /n=(1,16), settings
		endif
		settings[0][0]=SD_Dia
		settings[0][1]=SD_N
		settings[0][2]=SD_scat
		settings[0][3]=SD_Abs
		settings[0][4]=iter
		settings[0][5]=n_width
		settings[0][6]=k_width
		settings[0][7]=k_step
		settings[0][8]=n_step
		settings[0][10]=(n_width/n_step)*2*(k_width/k_step)*2
		settings[0][11]=spin_up_cycles
		settings[0][12]=spin_up_repeat
		settings[0][13]=abs_limit
		settings[0][14]=scat_limit
		settings[0][15]=Matrix_max



////runs code that estimates the probability of the discrete RI values test of having yielded the input retrieved RI.
		wave n_range, k_range
		RI_Interval (bin_input,k_input,N_input,Size_input,SD_Dia,SD_N,SD_Scat,SD_Abs,iter, wavelength,CL,calc_type)
		wave gauss_out


wave n_axis, k_axis

//calculates the number of idle/empty n and k bins relative to the number of bins
idle_n=gauss_out[mid_point[0]][31]/(numpnts(n_axis))
idle_k=gauss_out[mid_point[0]][32]/(numpnts(k_axis))

//prints the idle value
print "idle n ", gauss_out[mid_point[0]][31]/(numpnts(n_axis))
print "idle k ", gauss_out[mid_point[0]][32]/(numpnts(k_axis))






//measures the length of a single loop through the width adjustment and RI_interval function
		duration=datetime-start_time


		settings[0][9]=duration



		killwaves /z adaptive_output_all
		Make /n=(1,46), Adaptive_output_all


//records output from RI_interval in the adaptive_output_all matrix for each RI_interval run
		for (ad=0;ad<dimsize(gauss_out,1);ad+=1)
			Adaptive_output_all[0][ad]=gauss_out[mid_point[0]][ad]
		endfor





		if (adaptors<1)
			duplicate/o adaptive_output_all, output_tracker
		endif
		if (adaptors>0.1)
			concatenate/o/NP=0 {output_tracker,adaptive_output_all}, adaptive_output_all2
			duplicate /o adaptive_output_all2, output_tracker
		endif



		if (adaptors<1)

			duplicate/o settings, settings_tracker
		endif
		if (adaptors>0)
			concatenate/o/NP=0 {settings_tracker,settings}, settings_tracker2
			duplicate /o settings_tracker2, settings_Tracker
		endif




		variable tq



	endfor

setdatafolder ::


//creates an output wave including only the final solutions of RICE for each point in the input time series.

if (exists("adaptive_output")==0)
	Make /o /n=(dimsize(RI_t_series,0),46), Adaptive_output
endif
for (tq=0;tq<dimsize(gauss_out,1);tq+=1)
	Adaptive_output[iq][tq]=gauss_out[mid_point[0]][tq]
endfor




//creates simplified output wave and flag waves
	
	if(Adaptive_output[iq][26]<=1.05)
	space_n[iq]=0
	elseif (Adaptive_output[iq][26]<=1.125)
	space_n[iq]=1
	elseif (Adaptive_output[iq][26]<=1.2)
	space_n[iq]=2
	elseif (Adaptive_output[iq][26]>1.2)
	space_n[iq]=3
	endif



	if(Adaptive_output[iq][28]<=1.05)
	space_k[iq]=0
	elseif (Adaptive_output[iq][28]<=1.125)
	space_k[iq]=1
	elseif (Adaptive_output[iq][28]<=1.2)
	space_k[iq]=2
	elseif (Adaptive_output[iq][28]>1.2)
	space_k[iq]=3
	endif	


	if (Adaptive_output[iq][45]<=iter/2)
	count_[iq]=1
	else
	count_[iq]=0
	endif
	
	k_lower_bound[iq]=Adaptive_output[iq][24]
	k_upper_bound[iq]=Adaptive_output[iq][23]
	n_lower_bound[iq]=Adaptive_output[iq][22]
	n_upper_bound[iq]=Adaptive_output[iq][21]




	if (run_type==0)
	killdatafolder/Z $folder_name
	endif
	endif
	



endfor


killwaves /z RICE_Output	

duplicate /o Adaptive_output, RICE_Output

killwaves /z Adaptive_output

setdatafolder ::

end








///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//creates a confidence interval based on the sigmoidal fit of a probability distrubiton

Function Sigmoid2Interval(Sig_coef,step_width, CL)
wave Sig_coef
variable step_width,CL

variable base, Max_,center,Rate
Base=Sig_coef[0]
Max_=Sig_coef[1]	
Center=Sig_coef[2]
Rate=Sig_coef[3]


variable i=0, position=0, top, bottom,cuml_prob,x,y

variable non_neg_base
if (base<0)
	non_neg_base=0
	else
	non_neg_base=base
endif
i=0
for (i=0;cuml_prob<(1-((1-CL)/2));i+=1)


x=center+i*step_width
cuml_prob= (Base+(max_/(1+exp((center-x)/rate))))




top=x



if (i>100000)
cuml_prob=999
top=-999
endif


endfor




for (i=0;cuml_prob>((1-CL)/2);i+=1)
cuml_prob=0
x=0
x=center-i*step_width
cuml_prob=(Base+(max_/(1+exp((center-x)/rate))))

bottom=x


if (i>100000)
cuml_prob=-999
bottom=-999
endif

endfor

killwaves /z cI_out

make /n=(2) cI_out
if (bottom<0)
	bottom=0
endif

cI_out[0]=top
cI_out[1]=bottom
end


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Used to update the result of a previous k RICE run using a new confidence interval


Function Output_updater_k(Adaptive_output,CL)
wave adaptive_output
variable CL
make /o /n=(dimsize(adaptive_output,0)) k_lower_bound_updated
make /o /n=(dimsize(adaptive_output,0)) k_upper_bound_updated




variable i=0

make/o /n=(4,1) w_coef_1



for (i=0;i<dimsize(adaptive_output,0);i+=1)
	w_coef_1[0]=adaptive_output[i][41]
	w_coef_1[1]=adaptive_output[i][42]
	w_coef_1[2]=adaptive_output[i][43]
	w_coef_1[3]=adaptive_output[i][44]



	Sigmoid2Interval(w_coef_1,(w_coef_1[2]-0)/10000,CL)
	wave cI_out

k_upper_bound_updated[i]=cI_out[0]
k_lower_bound_updated[i]=cI_out[1]
endfor
killwaves w_coef_1
end

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Used to update the result of a previous n RICE run using a new confidence interval

Function Output_updater_n(Adaptive_output, CL)
wave adaptive_output
variable CL

make /o /n=(dimsize(adaptive_output,0)) n_lower_bound_updated
make /o /n=(dimsize(adaptive_output,0)) n_upper_bound_updated


variable i=0

make /o/n=(4,1) w_coef_1



for (i=0;i<dimsize(adaptive_output,0);i+=1)
	w_coef_1[0]=adaptive_output[i][37]
	w_coef_1[1]=adaptive_output[i][38]
	w_coef_1[2]=adaptive_output[i][39]
	w_coef_1[3]=adaptive_output[i][40]



	Sigmoid2Interval(w_coef_1,(w_coef_1[2]-1)/10000,CL)
	wave cI_out

n_upper_bound_updated[i]=cI_out[0]
n_lower_bound_updated[i]=cI_out[1]
endfor
killwaves w_coef_1, cI_out
end




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//makeRIrange creates an n and k input waves 
function makeRIrange(k_start,k_end,k_step,n_start,n_end,n_step)

variable k_start,k_end,k_step,n_start,n_end,n_step

variable n_points, k_points


if (exists ("matrix_max1")==1)
	wave matrix_max1
	n_points=(trunc(sqrt(matrix_max1[0])))+1
	k_points=(trunc(sqrt(matrix_max1[0])))+1
else


	n_points=abs((n_end-n_start)/(n_step))+2 //establishes the number of points needed for n and number needed of k
	k_points=abs((k_end-k_start)/(k_step))+1
endif


variable Rin=0, n=0,R
variable repeat=0, length=0

killwaves /z k_range, n_range

make  /o/n=(n_points), k_range    //creates n and k output waves
make /o /n=(n_points), n_range

duplicate/o n_range, n_holder
duplicate/o k_range, k_holder

for (r=0;r<numpnts(n_holder);r+=1)
	n_holder[r]=(n_step)*r+n_start
endfor


duplicate/o n_holder, n_range1
duplicate/o k_range, K_range1
K_range1=k_start

for (repeat=0;repeat<k_points;repeat+=1)
	k_holder=(k_step)*repeat+k_start

	if (repeat>0)
		concatenate /o/NP {n_range1,n_holder}, N_range
		concatenate /o/NP {k_range1,k_holder}, k_range
	
		duplicate/o k_range, K_range1
		duplicate/o n_range, n_range1
	endif

endfor

killwaves/z K_range1, n_range1, k_holder, n_holder

end




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//RI_Interval iterates inputs of a given size distribution and optical parameters using input uncertainties to estimates the distribution of RI values which could give a retrieved RI

function RI_Interval (Dia_bins,k_range,N_range,Size_profile,SD_Dia,SD_N,SD_Scat,SD_Abs,iter, wavelength,CL,calc_type) 
wave k_range,N_range,Size_profile, Dia_bins
variable SD_Dia,SD_N,SD_Scat,SD_Abs,iter, wavelength, CL,calc_type //establishes a threshold within which the result of an RI calcuations must be to the "obsereved" RI
variable abs_limit, scat_limit




//calculates the distance between tested n and k values
variable i=0,q=0,c=0, h=0, it=0,k_step,n_step

n_step=(n_range[1]-n_range[0])

for (i=1;i<numpnts(n_range);i+=1)
	if(n_range[i]-n_range[0]<(n_step/10))
		k_step=k_range[i]-k_range[0]
		i=numpnts(n_Range)
	endif
endfor
	
//creates a size distrubtion that has a point for each input RI value
killwaves /z size_dist
make /o /n=(numpnts(Dia_bins), numpnts(n_range)), size_dist

for (i=0;i<numpnts(Dia_bins);i+=1)
	for (c=0;c<numpnts(n_range);c+=1)
		size_dist[i][c]=size_profile[i]
	endfor
endfor


//Creates temporary bin diameter and number concentration waves to preserve the originals.	
duplicate /o dia_bins,Dia_input
duplicate/o size_profile, Num_input	



//creates a matrix for storage of the probability results
killwaves /z prob_dist
make /o /n=(numpnts(n_range),numpnts(n_range)) prob_dist
prob_dist=0
duplicate/o prob_dist, prob_holder

//calcuates the extinction, scattering, and absorption efficeny matrix for the base diameters
efficency_matrix(n_range, k_range, wavelength, dia_bins)  //establishes waves of the respective Qscat and Qabs for the size bins
wave Eff_matrix_k, Eff_matrix_n
duplicate /o Eff_Matrix_n, Eff_matrix_nbase
duplicate /o Eff_Matrix_k, Eff_matrix_kbase
conversion_2_dn(dia_bins)  // creates a dlog(d) wave for each bin
wave conversion
duplicate/o conversion, conversion_base  //saves the base dlog(d) matrix for the given bin_sizes
killwaves/z conversion
killwaves /z Eff_matrix_k, Eff_matrix_n


//calcuates the expected, unperturbed, optical coefficents 
Coef_Calc(dia_input,Eff_Matrix_nbase, Eff_matrix_kbase,size_profile) 
wave Ext_func_n, Ext_func_k


//calcuates the maximum observed difference between tested n and k values.  Sets the scattering and absorption reasonabilty thresholds to be within these values.  These thresholds are used to ensure neither scattering nor absorption are dominating the calculation.
Diff_calc_n_wave(Ext_func_n,n_Range) 
 wave max_obs_diff
 Scat_limit=max_obs_diff[0]*(0.5)
 
 if(scat_limit>2*(SD_Scat))
 	scat_limit=2*SD_scat
 endif
 
 
Diff_calc_k_wave(Ext_func_k,k_range)
 wave max_obs_diff
 abs_limit=max_obs_diff[0]*(0.5)
 
  if(abs_limit>2*(SD_abs))
 		abs_limit=2*SD_abs
 	endif
 


killwaves /z Ext_func_n,Ext_func_k
//pertubes the input values using the instrumental uncertainties.

for (it=0;it<iter;it+=1)
	variable bin_noise
	bin_noise=gnoise(SD_Dia)


dia_input=Dia_Bins*bin_noise+Dia_Bins


	for(q=0;q<numpnts(size_profile);q+=1)
		Num_input[q]=Size_profile[q]*gnoise(SD_N)+size_profile[q]
		if (num_input[q]<0)
			num_input[q]=0
		endif
	endfor


//calculates and efficency matrix using the pertubed bin locations
	efficency_matrix(n_range, k_range, wavelength, dia_input)  //establishes waves of the respective Qscat and Qabs for the size bins
	wave Eff_Matrix_n, Eff_matrix_k

	wave Eff_Matrix_n, Eff_matrix_k
	variable width=dimsize(size_dist,0)
	variable length=dimsize(size_dist,1)

	killwaves /z Sim_Scat, Sim_abs, conversion

	make /o /n=(numpnts(n_range)) Scat_input   //creates a matrix where rows will be the observed scattering at a given RI value for the corresponding size distrbution
	make /o /n=(numpnts(n_range)) Abs_Input   //creates a matrix where rows will be the observed absorption at a given RI value for the corresponding size distrbution

	
	//calcuates the optical coefficents using the perturbed diameters and number concentrations.
	Coef_Calc(dia_input,Eff_Matrix_n, Eff_matrix_k,Num_input)  //calculates the "tru" optical properties of the given size distrubtion over a given RI series
	
	//additionally perturbs optical coefficentsto account for observations uncertainties
	for (h=0;h<dimsize(n_range,0);h+=1)
		wave Ext_func_n, Ext_func_k
		Scat_input[h]=Ext_func_n[h]*gnoise(SD_Scat)+Ext_func_n[h]
		Abs_input[h]=Ext_func_k[h]*gnoise(SD_Abs)+Ext_func_k[h]
	endfor
	
	

	





	m_calc_choice(Scat_input, Dia_Bins, Wavelength, n_Range,k_range,Size_Dist,"results", abs_input,scat_limit,abs_limit,calc_type)


	wave RI_results, K_results


	variable as,ad,ag,af, n_in,k_in,prob_in

//assigns the results to a specific location given in the probabiliy distrubtion the calcated RI value (mr prime).  After many runs a probability distrubtion is created.
			
	for (as=0;as<numpnts(n_range);as+=1)
		n_in=RI_results[as]
		k_in=K_results[as]
		for(af=0;af<(numpnts(n_range));af+=1)
			if ((n_step/10)>(abs(n_in-n_range[af])))
				if ((k_step/10)>(abs(k_in-k_range[af])))
					Prob_Dist[as][af]= Prob_holder[as][af]+1
					duplicate/o prob_dist, prob_holder
				endif
			endif
		endfor
	endfor
	

	killwaves /z RI_results,K_results

endfor



MatrixOP /o Prob_sum=sumrows(prob_dist) 
duplicate /o prob_sum, mask_prob
mask_prob=1

//turns the probability distribtion from counts to decimal probabilities
for (i=0;i<dimsize(prob_dist,0);i+=1)
	for (c=0;c<dimsize(prob_dist,1);c+=1)
		prob_dist[i][c]=prob_holder[i][c]/iter
	endfor
endfor

wave mid_point


//analyzes the results using a sigmoidal function and calcuates a confidence interval around the input retrieved RI at the predefined confidence level.
n_k_dist_1d_simple(prob_dist,n_range,k_range,mid_point[0],prob_holder,CL)  
wave gauss_out


end








///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //creates a wave that will be used to convert the dn_dlogd matrix to a dn matrix

function conversion_2_dn(Dia_Bins)
	
wave Dia_Bins
duplicate /o Dia_bins, size_mean
 duplicate/O Size_mean, Conversion  
	//converts dn_dlog_d to dn
	
	variable i=0,size,q2=0, q_pre=0,size1,size2,size3
		variable q=0,p=0
		variable z=numpnts(Size_Mean)
		variable z2=z-1
			for(q=0;q<=z2;q+=1)     // loops through the midpoints and estimates the dlog(d) for a bin.  For the first and last bins, the distace between midpoints is assumed  
				q2=q+1
				q_pre=q-1
				size1=Size_mean[q]
				if (q==0)
					size2=Size_mean[q2]
					Conversion[q]=Log(size2)-Log(Size1)  //assumes the size of bin he first bin is the distance between midpoints of the first and second bin.
				else
					size3=Size_mean[q_pre]
					if (q<z2)
						size2=Size_mean[q2]
						Conversion[q]=(Log(Size2)-Log(Size1))/2+(Log(Size1)-Log(Size3))/2  // assumes the size of the middle bins are the sum of half the distance from the previous bin's midpoint and half the distance from the following bin's midpoint.
					endif
					if (q==z2)
						Conversion[q]=Log(Size1)-Log(Size3) //assumes the size of final bin is the distance between midpoints of the second to last bin and the last bin.
					endif
				endif	
			endfor
	
end



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


function Coef_Calc(bin_dia,Eff_Matrix_n, Eff_matrix_k,Size_Dist)
	wave bin_dia, Size_Dist,Eff_Matrix_n, Eff_matrix_k
	
	
	duplicate/o Eff_Matrix_n, output_abs
	duplicate/o Eff_Matrix_n, output_scat
	variable k=0
	variable g=0	
		variable i=0, ia=0
		variable q=0,p=0
		variable z=numpnts(bin_dia)
		variable z2=z-1
conversion_2_dn(bin_dia)
wave conversion
duplicate/o size_dist, dn_1
		
		
				dn_1=size_dist*conversion

	
	output_scat=0
		for(k=0;k<(dimsize(Eff_Matrix_n,1));k+=1)  // loops through the different n_vlaues
			for(g=0;g<(numpnts(bin_dia));g+=1) //loops through the different diameters
				output_scat[g][k]=Eff_Matrix_n[g][k]*dn_1(g)*((bin_dia(g)*(10^(-3)))^2)/4*pi  //  Calculates the extinctions for a given size bin(k) at a given time(g) for the given RI
			endfor
		endfor
	MatrixOP /O Ext_func_n = sumcols(output_scat)
	output_abs=0
			for(k=0;k<(dimsize(Eff_Matrix_k,1));k+=1)  // loops through the different n_vlaues
			for(g=0;g<(numpnts(bin_dia));g+=1) //loops through the different diameters
				output_abs[g][k]=Eff_Matrix_k[g][k]*dn_1(g)*((bin_dia(g)*(10^(-3)))^2)/4*pi  //  Calculates the extinctions for a given size bin(k) at a given time(g) for the given RI
			endfor
		endfor
	MatrixOP /O Ext_func_k = sumcols(output_abs)
end









 



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//preforms an optical closure calculation using the difference merit parameter



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//preforms an optical closure calculation using the chi squared merit parameter

function m_calc_choice(Scat, Size_Mean, Wavelength, RI_Range,k_range,count,results, abs_avg,scat_limit,abs_limit,calc_type)
	Wave Scat, Size_Mean, RI_Range,count, k_range, abs_avg
	variable Wavelength,scat_limit,abs_limit,calc_type
	string results
	duplicate/O Scat, Scat_avg //creates a scattering avg_wave, a relic from when had tried to include an averaging function in the macro
	variable numpts_RI=numpnts(RI_Range)	//establishes the length of the RI input value waves, for use in loops				
	variable numpts_size=numpnts(scat)  //establishes the length of the t_series, named numpts_size because this corresponds to the number of size distrubtion measurements
	
	
	make /O /N=(numpts_size) template  //creates waves where the output of the efficency calculations can be placed
	duplicate /O Size_Mean, Ext
	duplicate /O Size_Mean, Back
	duplicate /O count, output
	duplicate /O count, dn
	
	//converts dn_dlog_d to dn
	variable ia=0,t=0
	for (t=0;t<dimsize(count,0);t+=1)
		for (ia=0;ia<dimsize(count,1);ia+=1)	
			wave conversion_base
			dn[t][ia]=count[t][ia]*conversion_base[t]
		endfor
	endfor
	//dn=count
	//establishes waves for later functions
	duplicate /O count, output
	duplicate /O count, output_scat
	duplicate /O count, output_abs
	duplicate /O template,  Scat_series
	duplicate /O Size_Mean, Qe_
		
	output=0
	output_scat=0
	output_abs=0
	scat_series=0
		
	wave Eff_Matrix_nbase, Eff_matrix_kbase 
	//calculates the respective efficencies for each of the different diameters for the given RI
	make /o/n=(dimsize(size_mean,0)), abs_
	variable k=0
	variable g=0,nh=0

	variable dia_step, n_step, ext_n, abs_k, size_in,i



	for (i=0;i<numpnts(RI_range);i+=1)
		for(k=0;k<(numpnts(size_mean));k+=1)  // loops through the different diameters
			ext_n=Eff_Matrix_nbase[k][i]
			abs_k=Eff_Matrix_kbase[k][i]
			size_in=Size_mean[k]
			for(g=0;g<(numpts_size);g+=1) //loops through the different time poitns 
				output[k][g]=ext_n*dn[k][g]*((size_in*(10^(-3)))^2)/4*pi  //  Calculates the extinctions for a given size bin(k) at a given time(g) for the given RI
				output_abs[k][g]=abs_k*dn[k][g]*((Size_in*(10^(-3)))^2)/4*pi  //  Calculates the extinctions for a given size bin(k) at a given time(g) for the given RI
			endfor
		endfor
	
		
		MatrixOP /O Scat_Series = sumcols(output)   //sums the scattering at all sizebins for a given time point into a single number
		MatrixOP /O abs_Series = sumcols(output_abs)   //sums the scattering at all sizebins for a given time point into a single number
	
		// concatenate the expected scattering coeffiencts from each run of the for loop, eventually yeilding a RI vs time matrix with the values as the expected scattering
	
		variable f=i-1
		String name=num2str(i)+"_Summed_scat"  //names the sum
		String wavekiller="0_Summed_scat"
		String name_abs=num2str(i)+"_Summed_abs"
		String wavekiller2="0_Summed_abs"
	
		duplicate/O scat_series, $name  //renames scat_series by the number of the loop
		duplicate/O abs_series, $name_abs
		if (i==0)  //for the first one, don't concatenate just make new matrix
			duplicate /O $name,  theo_scat  //creates the wave theo_scat, which will have the final scattering expected for each RI input at each time point( a matrix of RI input vs time).
			duplicate /O $name_abs, theo_abs
	
		elseif(i>0)
			duplicate/O  theo_scat, disposable  //duplicates  theo_scat and creates a disposable wave so we can concatenate a new wave onto theo_Scat
			concatenate/O /NP=0{disposable,$name},  theo_scat  //adds the summed scattering to the theo_scat matrix
			killwaves $name  //cleans up unused waves
			duplicate/O theo_abs, disposable_abs
			concatenate/O /NP=0{disposable_abs,$name_abs}, theo_abs
			killwaves $name_abs 
		endif
		
		
	
	

		
	Endfor

	wave settings

	variable h=0
	//establish a series of waves that are used for the refractive index calculations
	duplicate/O scat_avg, RI_Results
	duplicate/O scat_avg, k_results
	duplicate/o scat_avg, Relative_Results
	duplicate/o scat_avg, Rel_abs
	duplicate/o scat_avg, Rel_scat
	
	k_results=-999  
	RI_results=-999
	variable uncr=1  // establishes and uncertainty variable, which will be used to estimate the range of RI's with a given scattering uncertainty(notice that this is built in, if you want to change it you have to change it here and at the bottom where the outputs are named)
	variable if1,if2,if3,if4,scat_diff, abs_diff,abs_obs,scat_obs,abs_delt,scat_delt,RI_Calc,Merit_Par,q
	//scattering difference calc
	q=dimsize(theo_scat,1)
	
	for (i=0; i<q;i+=1) //cycles through each timepoint
		Merit_Par=100
		abs_obs=abs_avg[i]
		scat_obs=scat[i]
		For (h=0; h<numpts_RI;h+=1) //cycles though each RI input point
			abs_diff=0
			scat_delt=abs(theo_scat[h][i]-scat_obs)
			scat_diff=abs((scat_delt/scat_obs))
			abs_delt=abs(Theo_abs[h][i]-abs_obs)
	
			if (abs_obs==0)
				abs_diff=abs_delt/0.000001
			else		
				abs_diff=abs((abs_delt)/abs_obs)
			endif
			
			if (calc_type==1)
			Merit_Par=((Scat_delt^2/(scat_obs*settings[0][2])^2)+(abs_delt^2/(abs_obs*settings[0][3])^2))
			elseif(calc_type==0)
			
			Merit_Par=Scat_delt+abs_delt
			endif
			
			if (h==0)
				RI_calc=	Merit_Par-0.00000001
			endif
				
			if (Merit_Par<RI_Calc)  //checks to see if the difference is the lowest one for that time point	
				rel_abs=abs_diff
				rel_scat=scat_diff
		 
				if(abs_diff<abs_limit)
					if(scat_diff<scat_limit)
						RI_Calc=Merit_Par
						 // if the difference is the lowest yet found, the dummy wave is set to that value
						RI_Results[i]=RI_Range[h]  // if the difference is the lowest yet found for timepoint , the n results is set to the corresponding n input value
						K_Results[i]=K_Range[h]  // if the difference is the lowest yet found timepoint, the k results is set to the corresponding k input value
					endif
				endif
				
			endif
				
				
		endfor
	endfor
	//DoWindow/K	name
	//display /n=name rel_abs, rel_scat
	// ModifyGraph rgb(Rel_abs)=(0,15872,65280)
	//abs difference calc, same as the scattering difference clac but for absorption
	//creates a dummy wave(RI_Calc) and results waves

			
		// the following script renames the by the input results string, and separtes results by uncertainty
	


			duplicate/O Relative_Results, Rel_results
	
	

	
	killwaves /z abs_Series, Back, combinedbarrow, disposable, disposable_abs, dn, ext, output, output_scat, output_abs, Qe_,Scat_avg, Scat_Series,template,theo_scat,theo_abs,merit_par,wheelbarrow_abs,absolutebarrow,absolutebarrow_abs 

	
	End



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Creates a diameter vs RI matrix of absorption efficiencies
function efficency_matrix(n_range, k_range wavelength, Dia_Bins) 
wave n_range, Dia_Bins,k_range


variable wavelength
variable i=0,q=0


duplicate/o dia_bins, Qe_
duplicate/o dia_bins, back
duplicate/o dia_bins, ext


make /O /N=(numpnts(Dia_Bins),numpnts(n_range)), Eff_Matrix_n
make /O /N=(numpnts(Dia_Bins),numpnts(n_range)), Eff_Matrix_k



for(i=0;i<numpnts(n_range);i+=1)
	variable RI=n_range(i), RI_k=k_range(i)
	Scatter_diameter(Dia_bins, cmplx(RI,RI_k),Wavelength, ext, Qe_, Back)
	
	for (q=0;q<numpnts(Dia_Bins);q+=1)
		Eff_Matrix_n[q][i]=Qe_[q]
		Eff_matrix_k[q][i]=Ext[q]-Qe_[q]
	endfor
endfor
end




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Creates the probablity distribution for n and k from the probability surface of m inputs vs m ouputs.  Fits the n and k probability distributions with sigmoids


function n_k_dist_1D_simple(prob_dist,n_range,k_range, mid_posit,prob_holder, CL)
wave prob_dist,n_range,k_range,prob_holder
variable mid_posit, CL
duplicate/o n_range, n_placement
duplicate/o k_range, k_placement

//newdatafolder/o n_k_prob_dists

variable i=0
for (i=0;i<numpnts(n_Range);i+=1)
	n_placement[i]=round((n_range[i]-n_range[0])/(n_range[1]-n_range[0]))
endfor

if (wavemax(k_range)>0)
	for (i=0;i<numpnts(n_Range);i+=1)
		k_placement[i]=round((k_range[i]-k_range[0])/(k_range[wavemax(n_placement)+1]-k_range[0]))
	endfor
else
	k_placement=0
endif
variable k_step=k_range[wavemax(n_placement)+1]-k_range[0]
variable n_step=n_range[1]-n_range[0]

make /o/n=(dimsize(prob_dist,1),46), gauss_out
make /o/n=(dimsize(prob_dist,1),2), fit_results

make /o/n=(round (wavemax(k_placement)+1)) k_dist
make /o/n=(round (wavemax(n_placement)+1)) n_dist

gauss_out=0
for (i=0;i<numpnts(n_dist);i+=1)
	n_dist[i]=n_range[i]
endfor
for (i=0;i<numpnts(k_dist);i+=1)
	k_dist[i]=k_range[i*round(wavemax(n_placement)+1)]
endfor

variable tot=0,totx=0,x,y




//for (tot=0;tot<dimsize(gauss_out,0);tot+=1)
tot=mid_posit
killstrings/z  strname
make /o/n=(round(wavemax(k_placement)+1),round(wavemax(n_placement)+1)) name, name2
name=0
for (totx=0;totx<numpnts(k_range);totx+=1)	
	
	
	x=k_placement[totx]
	y=n_placement[totx]
	
		
	name[x][y]=prob_dist[totx][tot]
	name2[x][y]=prob_holder[totx][tot]
	
endfor
duplicate/o name, name_holder
duplicate/o name, name_holder2
string strname="prob"+num2str(n_range[tot])+"_"+num2str(k_range[tot])
	
	

			
variable z=0
	
if	(sum(name_holder)>0)	
	name_holder2=name/(sum(name_holder))
else 
	name_holder2*=0
endif
gauss_out[tot][45]=sum(name2)
	
//	duplicate/o name_holder2, $strname

	
duplicate/o k_dist,k_axis
duplicate/o n_dist,n_axis
		
duplicate/o k_dist,k_prob
duplicate/o n_dist,n_prob	
matrixop /O n_prob1 = sumcols(name_holder2)
matrixop /O k_prob1 = sumrows(name_holder2)
variable	k_c, k_count=0,n_count=0
			
			
//transposes k_prob1 from columns to rows	
for(z=0;z<numpnts(k_dist);z+=1)
	k_prob[z]=k_prob1[z]
endfor

for(z=0;z<numpnts(n_dist);z+=1)
	n_prob[z]=n_prob1[0][z]
endfor

variable n_chi=0,k_chi=0, V_chisq=0, V_FitError=0, k_error =0, n_error=0

				
					
duplicate /o k_prob, k_avg
duplicate /o k_prob, k_dev
duplicate /o k_prob, k_cuml_dist

					
k_avg=k_prob*k_axis
variable k_sum=0,V_k_dev=0, k_median=0, k_density=0, k_empty_count=0, n_empty_count=0
k_sum=sum(k_avg)
k_dev=abs(k_axis-k_sum)*k_prob
V_k_dev=sum(k_dev)
variable	cum_dist_c=0, cnt=0, count_Stop
				
				
count_stop=0
for (cum_dist_c=0;cum_dist_c<numpnts(k_prob);cum_dist_c+=1)
				
	if( cum_dist_c>0)
		k_cuml_dist[cum_dist_c]=k_prob[cum_dist_c]+k_cuml_dist[cum_dist_c-1]
	else
		k_cuml_dist[cum_dist_c]=k_prob[cum_dist_c]
	endif
	if (k_median==0)
		if (k_cuml_dist[cum_dist_c]>0.5) 
			k_median=k_axis[cum_dist_c]
		endif
	endif
	if (k_density<k_prob[cum_dist_c])
		k_density=k_prob[cum_dist_c]				
		gauss_out[tot][30]=k_density
	ENDIF				
	if (count_Stop<1)
		if (k_cuml_dist[cum_dist_c]==0)
			k_empty_count+=1
		else
			count_stop=1			
		endif
	endif	
endfor
								
variable k_too_skinny=0
				
if (k_empty_count == 0)
	k_too_skinny=1
endif
				
count_stop=0
			
for (cum_dist_c=numpnts(k_prob)-1;cum_dist_c>-1;cum_dist_c-=1)	
	if (count_Stop<1)
		if (k_prob[cum_dist_c]== 0)
			k_empty_count+=1
			if (cum_dist_c==numpnts(k_prob)-1)
				k_too_skinny=1
			endif	
		else
		count_stop=1
		endif
	endif
endfor
				

gauss_out[tot][32]=k_empty_count
gauss_out[tot][34]=k_too_skinny
		
	
duplicate /o n_prob, n_avg
duplicate /o n_prob, n_dev
duplicate /o n_prob, n_cuml_dist
					
					
n_avg=n_prob*n_axis
variable n_sum=0,V_n_dev=0, n_median=0
n_sum=sum(n_avg)
n_dev=abs(n_axis-n_sum)*n_prob
V_n_dev=sum(n_dev)

cum_dist_c=0
cnt=0
variable n_denSIty=0
count_stop=0
for (cum_dist_c=0;cum_dist_c<numpnts(n_prob);cum_dist_c+=1)				
	if( cum_dist_c>0)
		n_cuml_dist[cum_dist_c]=n_prob[cum_dist_c]+n_cuml_dist[cum_dist_c-1]
	else
		n_cuml_dist[cum_dist_c]=n_prob[cum_dist_c]
	endif
					
	if (n_median==0)
		if (n_cuml_dist[cum_dist_c]>0.5) 
			n_median=n_axis[cum_dist_c]
		endif
	endif
					
	if (n_density<n_prob[cum_dist_c])
		n_density=n_prob[cum_dist_c]
		gauss_out[tot][29]=n_density 
	endif				
	if (n_cuml_dist[cum_dist_c]==0)
		n_empty_count+=1		
	else
		count_stop=1			
	endif					
endfor

variable n_too_skinny = 0
				
if (n_empty_count == 0)
	n_too_skinny=1
endif
count_stop=0
for (cum_dist_c=numpnts(n_prob)-1;cum_dist_c>-1;cum_dist_c-=1)	
	if (count_Stop<1)
		if (n_prob[cum_dist_c]==0)
			n_empty_count+=1
			if (cum_dist_c==numpnts(n_prob)-1)
				n_too_skinny=1
			endif
		else
			count_stop=1	
		endif
	endif	
endfor				
gauss_out[tot][31]=n_empty_count
gauss_out[tot][33]=n_too_skinny


				
string dist_name=num2str(n_range[tot])+"_"+num2str(k_range[tot])

print "number of k's", numpnts (k_axis)
print "number of n's", numpnts (n_axis)
	

//print "n high probability=", gauss_out[tot][29]
//print "k high probability=", gauss_out[tot][30]
	

				
wave final_run
					
if (final_run[0]==1)	
						
					
	killwaves /z w_sigma, w_coef
	V_FitError=0
	Variable V_fitoptions=4
	
	if(n_axis[0]-1<n_step)
		K0 = n_axis[0];K1 = 1;	
		CurveFit /H="1100"/M=2/W=0 /Q sigmoid, n_cuml_dist /X=n_axis /D		
	else
	K0 = 0;K1 = 1;	
	CurveFit /H="1100"/M=2/W=0 /Q sigmoid, n_cuml_dist /X=n_axis /D		
	endif
		
								
	if (exists("w_sigma")==1)
		wave w_sigma, w_coef
		duplicate/o w_coef, sig_n_coef_final
	else
		killwaves /z sig_n_coef_final
		make /o /n=4 sig_n_coef_final
		sig_n_coef_final=-999
	endif
							
									
	V_fitoptions=4
	V_FitError=0
	killwaves /z w_sigma, w_coef
	
	
	if(k_axis[0]-0<k_step)
	
	
	K0 = k_axis[0];K1 = 1;
	CurveFit /H="1100"/M=2/W=0 /Q sigmoid, k_cuml_dist /X=k_axis /D	
	else
	
	K0 = 0;K1 = 1;
	CurveFit /H="1100"/M=2/W=0 /Q sigmoid, k_cuml_dist /X=k_axis /D	
	endif



	if (exists("w_sigma")==1)
		wave w_sigma, w_coef	
 		duplicate/o w_coef, sig_k_coef_final 

	else
		killwaves /z sig_k_coef, sig_k_sigma_final
		make /o /n=4 sig_k_coef, sig_k_sigma_final
		duplicate/o sig_k_coef, sig_k_coef_final 
		sig_k_coef_final =-999
	endif
				
else
	make /o /n=4 sig_n_coef_final
	make /o /n=4 sig_k_coef_final
endif
							
							
killwaves /z w_sigma, w_coef
V_FitError=0
V_fitoptions=4	


if(n_axis[0]-1<n_step)


variable n_edge=1
K0 = n_axis[0];
killwaves /z w_sigma, w_coef
CurveFit /M=2/W=0 /Q sigmoid, n_cuml_dist /X=n_axis /D		

else
n_edge=0
killwaves /z w_sigma, w_coef
CurveFit /M=2/W=0 /Q sigmoid, n_cuml_dist /X=n_axis /D		
endif




								
if (exists("w_sigma")==1)
wave w_sigma, w_coef
duplicate/o w_coef, sig_n_coef
duplicate /o w_sigma, sig_n_sigma
else
killwaves /z sig_n_coef, sig_n_sigma
make /o /n=4 sig_n_coef, sig_n_sigma
sig_n_sigma=-999
sig_n_coef=-999
endif
							
V_fitoptions=4
V_FitError=0

if(k_axis[0]-1<k_step)
variable k_edge=1
K0 = k_axis[0];
killwaves /z w_sigma, w_coef
CurveFit /M=2/W=0 /Q sigmoid, k_cuml_dist /X=k_axis /D	
else
k_edge=0
killwaves /z w_sigma, w_coef
CurveFit /M=2/W=0 /Q sigmoid, k_cuml_dist /X=k_axis /D	
endif

if (exists("w_sigma")==1)
	wave w_sigma, w_coef	
	duplicate/o w_coef, sig_k_coef 
	duplicate /o w_sigma, sig_k_sigma
else
							
	duplicate/o sig_n_coef, sig_k_coef 
	duplicate /o sig_n_sigma, sig_k_sigma
							
	sig_k_coef =-999
	sig_k_sigma =-999
endif
				
							

gauss_out[tot][1]=sig_n_coef[0]
gauss_out[tot][2]=sig_n_sigma[0]
gauss_out[tot][3]=sig_n_coef[1]
gauss_out[tot][4]=sig_n_sigma[1]
gauss_out[tot][5]=sig_n_coef[2]
gauss_out[tot][6]=sig_n_sigma[2]
gauss_out[tot][7]=sig_n_coef[3]
gauss_out[tot][8]=sig_n_sigma[3]
gauss_out[tot][10]=sig_k_coef[0]
gauss_out[tot][11]=sig_k_sigma[0]
gauss_out[tot][12]=sig_k_coef[1]
gauss_out[tot][13]=sig_k_sigma[1]
gauss_out[tot][14]=sig_k_coef[2]
gauss_out[tot][15]=sig_k_sigma[2]
gauss_out[tot][16]=sig_k_coef[3]
gauss_out[tot][17]=sig_k_sigma[3]
gauss_out[tot][18]=V_n_dev
gauss_out[tot][19]=n_median
gauss_out[tot][20]=k_median
gauss_out[tot][21]=V_k_dev
gauss_out[tot][37]=sig_n_coef_final[0]
gauss_out[tot][38]=sig_n_coef_final[1]
gauss_out[tot][39]=sig_n_coef_final[2]
gauss_out[tot][40]=sig_n_coef_final[3]	
gauss_out[tot][41]=sig_k_coef_final[0]
gauss_out[tot][42]=sig_k_coef_final[1]
gauss_out[tot][43]=sig_k_coef_final[2]
gauss_out[tot][44]=sig_k_coef_final[3]	
		
		
if (gauss_out[tot][1]+(1-gauss_out[tot][3])>0.05)
	gauss_out[tot][0]=1
endif	
				
if (final_run[0]==0)
	Sigmoid2Interval(sig_n_coef,n_step/1000,CL)			
else
	Sigmoid2Interval(sig_n_coef_final,n_step/1000,CL)
endif
wave cI_out
gauss_out[tot][21]=cI_out[0]
gauss_out[tot][22]=cI_out[1]		

if (cI_out[1]<1)	
	gauss_out[tot][22]=1
endif

if (final_run[0]==0)
	Sigmoid2Interval(sig_k_coef,K_step/1000,CL)
else
	Sigmoid2Interval(sig_k_coef_final,K_step/1000,CL)
endif
	wave cI_out
	gauss_out[tot][23]=cI_out[0]
	gauss_out[tot][24]=cI_out[1]
	
	if (cI_out[1]<0)
	
		gauss_out[tot][24]=0
	endif
							
if ((sig_n_sigma[2]/sig_n_coef[2])>0.2)
	gauss_out[tot][25]=3
elseif ((sig_n_sigma[2]/sig_n_coef[2])>0.01)
	gauss_out[tot][25]=2
elseif ((sig_n_sigma[2]/sig_n_coef[2])>0.001)
				
	gauss_out[tot][25]=1
else
	gauss_out[tot][25]=0
endif				
				
//print "n fit", gauss_out[tot][25]


if(n_edge==0)
						
	if (sig_n_coef[1]-sig_n_coef[0]>1.2)
		gauss_out[tot][26]=3
	elseif (sig_n_coef[1]-sig_n_coef[0]>1.125)
		gauss_out[tot][26]=2
	elseif (sig_n_coef[1]-sig_n_coef[0]>1.05)
		gauss_out[tot][26]=1
	else
		gauss_out[tot][26]=0
	endif				
else
	if (sig_n_coef[1]-sig_n_coef[0]+n_prob[0]>1.2)
		gauss_out[tot][26]=3
	elseif (sig_n_coef[1]-sig_n_coef[0]+n_prob[0]>1.125)
		gauss_out[tot][26]=2
	elseif (sig_n_coef[1]-sig_n_coef[0]+n_prob[0]>1.05)
		gauss_out[tot][26]=1
	else
		gauss_out[tot][26]=0
	endif
endif


			
print "n Height",gauss_out[tot][26]




	
	
	

if (sig_k_sigma[2]/sig_k_coef[2]>0.2)
	gauss_out[tot][27]=3
elseif (sig_k_sigma[2]/sig_k_coef[2]>0.01)
	gauss_out[tot][27]=2
elseif (sig_k_sigma[2]/sig_k_coef[2]>0.005)
	gauss_out[tot][27]=1
else
	gauss_out[tot][27]=0
endif				
			
//print "k fit=", gauss_out[tot][27]


if (k_edge==0)
	if (sig_k_coef[1]-sig_k_coef[0]>1.2)
				
		gauss_out[tot][28]=3
	elseif (sig_k_coef[1]-sig_k_coef[0]>1.125)
		gauss_out[tot][28]=2
	elseif (sig_k_coef[1]-sig_k_coef[0]>1.05)
		gauss_out[tot][28]=1
	else
		gauss_out[tot][28]=0
	endif		

else	
	if (sig_k_coef[1]-sig_k_coef[0]+k_prob[0]>1.2)
				
		gauss_out[tot][28]=3
	elseif (sig_k_coef[1]-sig_k_coef[0]+k_prob[0]>1.125)
		gauss_out[tot][28]=2
	elseif (sig_k_coef[1]-sig_k_coef[0]+k_prob[0]>1.05)
		gauss_out[tot][28]=1
	else
		gauss_out[tot][28]=0
	endif	
endif
				
print "k Height", gauss_out[tot][28]
				
				
				
				 
									 
duplicate /o n_prob, n_real
duplicate/o k_prob, k_real
duplicate /o n_cuml_dist, n_real_cdist
duplicate /o n_cuml_dist, n_real_cdist_high
duplicate /o n_cuml_dist, n_real_cdist_low
duplicate /o k_cuml_dist, k_real_cdist
duplicate /o k_cuml_dist, k_real_cdist_high
duplicate /o k_cuml_dist, k_real_cdist_low
killwaves /z n_name_holder, name_holder_2


MatrixOP /o Prob_sum=sumrows(prob_holder) 
 gauss_out[tot][35]=prob_sum[tot]


		
end		










///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Calculates the difference between scattering coefficents for the tested refractive index values

Function Diff_calc_n_wave(Coef_wave,RI_wave)
wave Coef_wave, RI_wave
killwaves/z max_obs_diff, coef_input
make /n=1 max_obs_diff
variable pnt_cnt=0,obs_diff=0 


if (dimsize(coef_wave,1)<dimsize(coef_wave,0))
	matrixop coef_input=coef_Wave^t
else
duplicate /o coef_wave, coef_input
endif


for(pnt_cnt=1;pnt_cnt<dimsize(coef_input,1);pnt_cnt+=1)

if(RI_wave[pnt_cnt]>RI_wave[pnt_cnt-1])
if(coef_input[pnt_cnt-1]>0)
	obs_diff=(coef_input[0][pnt_cnt]-coef_input[0][pnt_cnt-1])/coef_input[0][pnt_cnt-1]
	if (obs_diff>(max_obs_diff[0]))
		max_obs_diff[0]=obs_diff
	endif
endif
endif
endfor

killwaves coef_input
	
end




	


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Calculates the difference between absorption coefficents for the tested refractive index values
Function Diff_calc_k_wave(Coef_wave,RI_wave)
wave Coef_wave, RI_wave
killwaves/z max_obs_diff, coef_input
make /n=1 max_obs_diff
variable pnt_cnt=0,obs_diff=0, interval_size
 
 if (dimsize(coef_wave,1)<dimsize(coef_wave,0))
	matrixop coef_input=coef_Wave^t
else
duplicate /o coef_wave, coef_input
endif
 
 

for(pnt_cnt=1;pnt_cnt<dimsize(coef_input,1);pnt_cnt+=1)



if(RI_wave[pnt_cnt]-RI_wave[pnt_cnt-1]>0)

interval_size=pnt_cnt
pnt_cnt=dimsize(coef_input,1)
endif
endfor


for(pnt_cnt=interval_size;pnt_cnt<dimsize(coef_input,1);pnt_cnt+=1)
if(coef_input[pnt_cnt-1]>0)
if ((coef_input[0][pnt_cnt-interval_size])>0)
	obs_diff=(coef_input[0][pnt_cnt]-coef_input[0][pnt_cnt-interval_size])/coef_input[0][pnt_cnt-interval_size]
	
	if (obs_diff>(max_obs_diff[0]))
		max_obs_diff[0]=obs_diff
	endif
endif
endif
endfor
	
killwaves coef_input	
end

  

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//bin_trim removes bin's within the size distrubtion that contribute less than 0.01% of the total absorption and scattering
function bin_trim(n_value,k_Value,n_width,k_width,dia_bins,size_profile1,wavelength)
wave dia_bins,size_profile1
variable n_value,k_value,n_width,k_width,wavelength

variable ty,tx,bin_counter,bin_counter2,qi

makeRIrange(k_value-k_width,k_value+k_width,k_width,n_value-n_width,n_value+n_width,n_width)  //makes a RI ranges

wave n_Range, k_range
efficency_matrix(n_range, k_range, wavelength, dia_bins)  //establishes waves of the respective Qscat and Qabs for the size bins

wave Eff_matrix_k, Eff_matrix_n


Coef_Calc(dia_bins,Eff_Matrix_n, Eff_matrix_k,size_profile1) //calcuates the total Optical coefficents and the portiton of the optical coefficents associated with each bin at each RI tested
killwaves /z Eff_matrix_k, Eff_matrix_n
wave Ext_func_n, Ext_func_k, output_abs, output_scat


for(ty=0;ty<dimsize(output_abs,1);ty+=1)
	for(tx=0;tx<dimsize(output_abs,0);tx+=1)
		output_abs[tx][ty]/=ext_func_k[ty]
		output_scat[tx][ty]/=ext_func_n[ty]
	endfor
endfor


MatrixOP /o output_abs_t=output_abs^t
MatrixOP /o output_scat_t=output_scat^t



MatrixOP /o output_abs_t=output_abs^t
MatrixOP /o output_scat_t=output_scat^t


make /o/n=(dimsize(output_abs_t,1)) max_abs_bin, max_scat_bin
variable opmax1, opmax2
for(opmax1=0;opmax1<dimsize(output_abs_t,1);opmax1+=1)
	for(opmax2=0;opmax2<dimsize(output_abs_t,0);opmax2+=1)
	
	if (output_abs_t[opmax2][opmax1]>max_abs_bin[opmax1])
max_abs_bin[opmax1]=output_abs_t[opmax2][opmax1]
	endif
	if (output_scat_t[opmax2][opmax1]>max_scat_bin[opmax1])
max_scat_bin[opmax1]=output_scat_t[opmax2][opmax1]
	endif
endfor
endfor


/////////////Creates a size profile which only uses potentiall important diameter bins so the whole analysis goes faster



bin_counter=0
bin_counter2=-1

duplicate /o size_profile1, size_profile2, bin_input
size_profile2=0
bin_input=0

for(qi=0;qi<numpnts(dia_bins);qi+=1)
	if	(max_abs_bin[qi][0]>0.01)
		
		bin_counter=bin_counter2+1
		bin_counter2=bin_counter
		size_profile2[bin_counter]=size_profile1[qi]
		bin_input[bin_counter]=dia_bins[qi]
		
	elseif	(max_scat_bin[qi][0]>0.01)
	
	
		bin_counter=bin_counter2+1
		bin_counter2=bin_counter
		size_profile2[bin_counter]=size_profile1[qi]
		bin_input[bin_counter]=dia_bins[qi]
		
	endif		
endfor


redimension /N=(bin_counter+1), size_profile2
redimension /N=(bin_counter+1), bin_input
 
duplicate /o size_profile2, size_input



killwaves /z size_profile1, output_abs,output_scat,max_Abs_bin,max_scat_bin,Ext_func_n, Ext_func_k, n_Range, k_range, Ext_func_n, Ext_func_k
 end
 
 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 //Originally Created by Charles A Brock, NOAA- ESRL Chemical Sciences Division, Boulder, CO 80305
 //     Calls BHMIE for a set of sizes 
//	Requires  diam, and refractive index wave and wavelength.
//  	as function of diam.  Returns waves containing Qext,Qscat, and Qback as f(diameter),
//	which can then be integrated across the size distribution.
//
      function Scatter_diameter(diam,ref_ind,wlngth,ext,scat,back)
        wave diam,ext,scat,back
        variable wlngth
        variable /c ref_ind
                
        variable ndiam=numpnts(diam),j,x,nang=5
        make/o/n=3 QOpt
        wave QOpt
                
         j=-1
        do
        	j+=1
		X= pi*diam[j]/wlngth
	//	ref_ind=ref_ind_wave[j]    
	//	print j,X,ref_ind,nang  	
        	BHMie (X, ref_ind, nang, QOpt)
        	//print QOpt[0],QOpt[1],Qopt[2]
        	scat[j]=QOpt[0]
        	ext[j]=QOpt[1]
        	back[j]=QOpt[2]
        while (j<ndiam-1)
        
      
       Return 0
	end
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
    /BHMie
//
//
//     SUBROUTINE BHMIE CALCULATES AMPLITUDE SCATTERING MATRIX ELEMENTS
//     AND EFFICIENCIES FOR EXTINCTION, TOTAL SCATTERING, AND BACKSCAT-
//     TERING FOR A GIVEN SIZE PARAMETER  AND RELATIVE REFRACTIVE INDEX
//
      function BHMie (X, REFREL, NANG, QOPT)
      variable x,nang
      variable /d/c refrel
      wave QOpt
      
      make /o/d/n=300 amu,theta,tau,PI0,PI1,PII
      make /d/o/c/n=3000 D=cmplx(0,0)
      make /d/o/c/n=(600) S1=cmplx(0,0),S2=cmplx(0,0)
      variable /d/g/c Y=cmplx(0,0)
      variable /d/g/c XI=cmplx(0,0),XI0=cmplx(0,0),XI1=cmplx(0,0),AN=cmplx(0,0),BN=cmplx(0,0)
      variable /d DX,Xstop,Ymod,dang,PSI0,PSI1,CHI0,CHI1,APSI0, APSI1,PP,TT,PSI,APSI,CHI
      variable  I,J,K,N,NN,RN,DN,FN,JJ,Nstop,NMX
	variable /d qext,qsca,qback

      DX = X
      Y = X * REFREL
//
//     SERIES TERMINATED AFTER NSTOP TERMS
//
      XSTOP = X + 4. * X ^ 0.3333 + 2.0
      NSTOP = floor(XSTOP)
      YMOD = CABS(Y)
      NMX =ceil(MAX(XSTOP,YMOD)) + 15
      DANG = (pi/2.)/ (NANG - 1)
      J=0
      DO //555 J = 1, NANG
        J+=1
        THETA[J] = (J-1.) * DANG
	   AMU[J] = COS (THETA[J])
 	While (J<Nang)   //555

//     LOGARITHMIC DERIVATIVE D(J) CALCULATED BY DOWNWARD RECURRENCE
//     BEGINNING WITH INITIAL VALUE 0.0 + I*0.0 AT J = NMX
//
      D[NMX] = CMPLX (0.0,0.0)
      NN = NMX - 1

      N=0
      DO //120
        N+=1
        RN = NMX - N + 1
	  D[NMX - N] = (RN/Y) - (1./(D[RN] + RN/Y))
	  // 120
	While (N<NN)  

	J=0
      DO   //666
	  J+=1
        PI0[J] = 0.0
        PI1[J] = 1.0
        // 666
      while(J<NANG)
      
      NN = 2 * NANG - 1
      J=0
      DO  //777
        J+=1
        S1[J] = CMPLX(0.0,0.0)
        S2[J] = CMPLX(0.0,0.0)
      //777
      While (J<NN)
 	
//
//     RICATTI-BESSEL FUNCTIONS WITH REAL ARGUMENT X CALCULATED BY
//     UPWARD RECURRENCE
//
      PSI0 = COS (DX)
      PSI1 = SIN (DX)
      CHI0 = -1*SIN (X)
      CHI1 = COS (X)
      APSI0 = PSI0
      APSI1 = PSI1
      XI0 = CMPLX(APSI0,-1*CHI0)
      XI1 = CMPLX(APSI1,-1*CHI1)
      QSCA = 0.0
      N = 1
      
      DO  // 200
        DN = N
        RN = N
        FN = (2. * RN + 1.) / (RN * (RN + 1.))
        PSI = (2. * DN - 1.) * PSI1 / DX - PSI0
        APSI = PSI
        CHI = (2. * RN - 1.) * CHI1 / X - CHI0
        XI = CMPLX (APSI, -1*CHI)
        AN = (D[N] / REFREL + RN / X) * APSI - APSI1
        AN = AN / ((D[N] / REFREL + RN/X) * XI - XI1)
        BN = (REFREL * D[N] + RN/X) * APSI - APSI1
        BN = BN / ((REFREL * D[N] + RN/X) * XI - XI1)
        QSCA=QSCA+(2.*RN+1.)*(CABS(AN)*CABS(AN)+CABS(BN)*CABS(BN))
  
        J=0
        DO //789 J = 1, NANG
           J+=1
           JJ = 2 * NANG - J
           PII[J] = PI1[J]
           TAU[J] = RN * AMU[J] * PII[J] - (RN + 1.) * PI0[J]
           PP = (-1.)^(N-1)
           S1[J] = S1[J] + FN * (AN * PII[J] + BN * TAU[J])
           TT= (-1.)^N
           S2[J] = S2[J] + FN * (AN * TAU[J] + BN * PII[J])
           IF (J  != JJ)
             S1[JJ] = S1[JJ] + FN * (AN * PII[J] * PP + BN * TAU[J] * TT)
             S2[JJ] = S2[JJ] + FN * (AN * TAU[J] * TT + BN * PII[J] * PP)
           EndIf
	    // 789  
 	  While (J<Nang)
	
        PSI0 = PSI1
        PSI1 = PSI
        APSI1 = PSI1
        CHI0 = CHI1
        CHI1 = CHI
        XI1 = CMPLX(APSI1, -1*CHI1)
        N = N + 1
        RN = N
        J=0
        DO //999 J = 1, NANG
          J+=1
          PI1[J] = ((2. * RN - 1.) / (RN - 1.)) * AMU[J] * PII[J]
          PI1[J] = PI1[J] - RN * PI0[J] / (RN - 1.)
          PI0[J] = PII[J]
          //999
        While (J<Nang)
 
      While ((N-1-NSTOP) < 0) 

	QSCA = (2. / (X*X)) * QSCA
	QEXT = (4. / (X*X)) * REAL (S1[1])
	QBACK = (4./ (X*X)) * CABS(S1[2 * NANG - 1]) * CABS(S1[2*NANG-1])
	QOpt[0]=QSca
	QOpt[1]=QExt
	QOpt[2]=QBack
//
//
	//print QSCA,QEXT,QBACK
	RETURN qext
	End Function
	
	
	
	
	
	function RI_calculation()  //actual macro
	
	String nScat="avg_Scat", Size_Mean="bin_dia", Count="size_dist", Results,nabs_avg="avg_abs",Scat_SD="Scat_SD",Abs_sd="Abs_SD"
	

variable k_step=0.001, k_end=0.02, k_start=0, n_start=1, n_end=1.8, n_step=0.01, Thres=0.1, T_displace=0 ,min_diff_scat=0.1,Min_diff_ABS=0.1,calc_type=0

prompt k_start, "lowest k value to test"
prompt k_end, "highest k value to test"
prompt k_step, "width of k steps"
prompt n_start, "Lowest n value to test"
prompt n_end, "highest n value to test"
prompt n_step, "width of n steps"
prompt Scat_SD, "Standard Deviation Wave of Scattering Measurements", popup, WaveList("*",";","")
prompt Abs_SD,  "Standard Deviation Wave of Absorption Measurements", popup, WaveList("*",";","")
Prompt calc_type, "RI Calcuation type, 0=difference based (Frie et al. 2019), 1=\"chi squared based\" (Dingle et al. 2019)"
doprompt "Establish RI interation range", k_start,k_end,k_step,n_start,n_end,n_step


doprompt "Calcuation Details" scat_sd,abs_sd,calc_type


if(v_flag==1) //aborts on cancel
abort 
print "canceled by user"
endif


makeRIrange(k_start,k_end,k_step,n_start,n_end,n_step)
wave n_range, k_range
wave RI_range=n_range






	variable Wavelength=375,uncertainty,numpts_RI,DL_scat=1.39, DL_Abs=1.0
	prompt nScat, "Observed Scattering Time Series",popup, WaveList("*",";","")
	prompt nabs_Avg "Observed Absorption Time Series" ,popup, WaveList("*",";","")
	Prompt DL_scat, "Scattering Detection Limit"
	Prompt DL_Abs, "Absorption detection limit"
	Prompt Thres, "Maximum relative acceptable change in mode between consecutive size measurements" 
	prompt Size_Mean, "Midpoint diameter's of size bins (nm)?",popup, WaveList("*",";","")
	prompt Count, "dn dlog(d) matrix",popup, WaveList("*",";","")
	prompt Wavelength, "Wavelength of Optical Measurements"   
	doprompt "Data Input", nScat, nabs_avg,size_mean,count,wavelength,DL_scat,DL_abs,Thres
	
if(v_flag==1) //aborts on cancel
abort 
print "canceled by user"
endif
	


numpts_RI=numpnts(RI_Range)
 
 
wave abs_avg=$nabs_Avg, scat=$nscat
variable i
 
for (i=0;i<numpnts(scat);i+=1)  //removes BDL data from averages, DLs can be set above
	if (scat[i]<DL_scat)
		scat[i]=nan
	endif
	if (abs_Avg[i]<DL_abs)
		abs_Avg[i]=0
	endif
endfor
 
 
 

m_calc_simple(Scat,abs_avg, $Size_Mean, Wavelength, RI_Range,K_Range,$count,$Scat_SD,$Abs_SD,Calc_type)  // Calls function to calculate the RI

wave results_n,results_k

//DetermineMode($count, $size_mean)
//wave Max_Diam

//RI_filter(thres,results_n, results_k, max_diam)

//wave results_k_mask, results_n_mask
//duplicate /o results_k_mask, results_k
//duplicate /o results_n_mask, results_n

edit results_k,results_n


killwaves /z amu,D,PI1,PI0,QOpt,S1,S2,tau,theta,PII,RI_results, results_k_mask, results_n_mask  // cleans up waves created by functions
End




function m_calc_simple(Scat,abs_avg, Size_Mean, Wavelength, RI_Range,k_range,count,Scat_SD,Abs_SD,Calc_type)
	Wave Scat, Size_Mean, RI_Range,count, k_range, abs_avg,Scat_SD,Abs_SD
	variable Wavelength,Calc_type
	

	
	duplicate/O Scat, Scat_avg //creates a scattering avg_wave, a relic from when had tried to include an averaging function in the macro
	variable numpts_RI=numpnts(RI_Range)	//establishes the length of the RI input value waves, for use in loops				
	variable numpts_size=numpnts(scat)  //establishes the length of the t_series, named numpts_size because this corresponds to the number of size distrubtion measurements
	
	
	make /O /N=(numpts_size) template  //creates waves where the output of the efficency calculations can be placed
	duplicate /O Size_Mean, Ext
	duplicate /O Size_Mean, Back
	duplicate /O count, output
	duplicate /O count, dn
	duplicate/O Size_mean, Conversion  //creates a wave that will be used to convert the dn_dlogd matrix to a dn matrix
	
	//converts dn_dlog_d to dn
	variable i=0
		variable q=0,p=0
		variable z=numpnts(Size_Mean)
		variable z2=z-1
			for(q=0;q<=z2;q+=1)     // loops through the midpoints and estimates the dlog(d) for a bin.  For the first and last bins, the distace between midpoints is assumed  
				variable q2=q+1, q_pre=q-1
				if (q==0)
					Conversion[q]=Log(Size_mean(q2))-Log(Size_mean(q))  //assumes the size of bins is the distance between midpoints for the first bin
					else
						if (q<z2)
							Conversion[q]=(Log(Size_mean(q2))-Log(Size_mean(q)))/2+(Log(Size_mean(q))-Log(Size_mean(q_pre)))/2  // assumes the size of the size bin is the sume of half the distance from the previous bin and half from the next bin.
						endif
						if (q==z2)
							Conversion[q]=Log(Size_mean(q))-Log(Size_mean(q_pre)) //assumes the size of bins is the distance between midpoints for the last bin
						endif
				endif	
			endfor
			for(q=0;q<(numpts_size);q+=1)  // creates the dn matrix by multiplying the conversion factor by the dn_dlog(d) matrix
				for(p	=0;p<(z);p+=1)
					dn[p][q] =Conversion(p)*count[p][q]
				endfor
			endfor
		//establishes waves for later functions
		duplicate /O count, output
		duplicate /O count, output_scat
		duplicate /O count, output_abs
		duplicate /O template,  Scat_series
		duplicate /O Size_Mean, Qe_
		
	for (i=0;i<numpts_RI;i+=1)  // starts for loop for efficency calcs and establishes the input n and k for said calc
		variable x=RI_Range[i]
		variable RI_k=k_Range[i]
	
	Scatter_diameter(Size_mean, cmplx(x,RI_k),Wavelength, Ext, Qe_, Back) //calculates the respective efficencies for each of the different diameters for the given RI
	variable k=0
	variable g=0
		for(k=0;k<(z);k+=1)  // loops through the different diameters
			for(g=0;g<(numpts_size);g+=1) //loops through the different time poitns 
				output[k][g]=Ext(k)*dn[k][g]*((Size_mean(k)*(10^(-3)))^2)/4*pi  //  Calculates the extinctions for a given size bin(k) at a given time(g) for the given RI
				output_scat[k][g]=Qe_(k)*dn[k][g]*((Size_mean(k)*(10^(-3)))^2)/4*pi  //  Calculates the scattering for a given size bin(k) at a given time(g) for the given RI
				output_abs[k][g]=output[k][g]-output_scat[k][g]  //  Calculates the absorption for a given size bin(k) at a given time(g) for the given RI
			endfor
		endfor
		
		
	MatrixOP /O Scat_Series = sumcols(output_scat)   //sums the scattering at all sizebins for a given time point into a single number
	
	// concatenate the expected scattering coeffiencts from each run of the for loop, eventually yeilding a RI vs time matrix with the values as the expected scattering
	
	variable f=i-1
	String name=num2str(i)+"_Summed_scat"  //names the sum
	String wavekiller="0_Summed_scat"
	
	duplicate/O scat_series, $name  //renames scat_series by the number of the loop
	
		if (i==0)  //for the first one, don't concatenate just make new matrix
			duplicate /O $name, Trowel  //creates the wave trowel, which will have the final scattering expected for each RI input at each time point( a matrix of RI input vs time).
		endif
		if(i>0)
			duplicate/O Trowel, disposable  //duplicates the trowel and creates a disposable wave so we can concatenate a new wave onto trowel (
			concatenate/O /NP=0{disposable,$name}, Trowel  //adds the summed scattering to the trowel matrix
			killwaves $name  //cleans up unused waves
		endif
		
		
		MatrixOP /O abs_Series = sumcols(output_abs)  //same as above, except for absorption
	//scatering concatenate the scattering coefficents tseries 
	String name_abs=num2str(i)+"_Summed_abs"
	String wavekiller2="0_Summed_abs"
	
	duplicate/O abs_series, $name_abs
	
		if (i==0)
			duplicate /O $name_abs, Trowel_abs
		endif
		if(i>0)
			duplicate/O Trowel_abs, disposable_abs
			concatenate/O /NP=0{disposable_abs,$name_abs}, Trowel_abs
			killwaves $name_abs 
		endif
		
Endfor
killwaves $wavekiller,$wavekiller2
	variable h=0
	//establish a series of "wheelbarrows" (get it? they carry things!) that are used for the refractive index calculations
	duplicate/o Trowel, WheelBarrow
	duplicate/o Trowel, AbsoluteBarrow
	duplicate/o Trowel, WheelBarrow_Abs
	duplicate/o Trowel_abs, AbsoluteBarrow_abs
	duplicate/o Trowel, combinedbarrow
	
	
	//scattering difference calc
	
		If (calc_type==0)
		for (i=0; i<q;i+=1)
			For (h=0; h<numpts_RI;h+=1)
				WheelBarrow[h][i]=Trowel[h][i]-Scat_avg[i] //calculates the difference between the calculate Bscat and the observed Bscat and puts in the the wheelbarrow
			endfor
		endfor
	
		elseif (calc_type==1)
				for (i=0; i<q;i+=1)
					For (h=0; h<numpts_RI;h+=1)
						WheelBarrow[h][i]=(Trowel[h][i]-scat_avg[i])^2/(scat_SD[i]^2)
					endfor
				endfor
		endif
		
		for (i=0; i<q;i+=1)
			For (h=0; h<numpts_RI;h+=1)
				AbsoluteBarrow[h][i]=((WheelBarrow[h][i])^2)^(1/2)  //takes the wheelbarrow and finds the absolute values of the differences and puts it in the absolutebarrow
			Endfor
		Endfor
	
	//abs difference calc, same as the scattering difference clac but for absorption
		If (calc_type==0)
			for (i=0; i<q;i+=1)
				For (h=0; h<numpts_RI;h+=1)
					WheelBarrow_abs[h][i]=Trowel_abs[h][i]-abs_avg[i]
				endfor
			endfor
		
		elseif (calc_type==1)
				for (i=0; i<q;i+=1)
					For (h=0; h<numpts_RI;h+=1)
						WheelBarrow_abs[h][i]=(Trowel_abs[h][i]-abs_avg[i])^2/(abs_SD[i]^2)
					endfor
				endfor
		endif
		
		
	
		for (i=0; i<q;i+=1)
			For (h=0; h<numpts_RI;h+=1)
				AbsoluteBarrow_Abs[h][i]=((WheelBarrow_abs[h][i])^2)^(1/2)  // same as above but for a  absorption
			Endfor
		Endfor
	
	//creates a dummy wave(RI_Calc) and results waves
	duplicate/O Scat_Avg, RI_Calc
	duplicate/O scat_avg, RI_Results
	duplicate/O scat_avg, k_results
	
	k_results=scat_avg*0  //sets output values to 0, so if there are nan values they become 0s
	RI_results=scat_avg*0
	
		RI_Calc=Scat_Avg+10000000000	// sets the dummy wave to a very high value
		
		
	
		
		CombinedBarrow=absolutebarrow+absolutebarrow_Abs  //adds the absolute difference of the scattering and absorption together, to get a summed difference acounting for both scattering and abs
		

		variable pnt,scat_limit, abs_limit
		killwaves /z possible_abs,possible_scat

		make /n=(dimsize(absolutebarrow,0)) possible_abs,possible_scat
		for (i=0; i<q;i+=1) //cycles through each timepoint
	
			 
			 for (pnt=1;pnt<(dimsize(absolutebarrow,0));pnt+=1)
			 		possible_abs[pnt]=Trowel_abs[pnt][i]
			 		possible_scat[pnt]=Trowel[pnt][i]
			 endfor
			 
			 
			  Diff_calc_n_wave(possible_scat,RI_range) 
			wave max_obs_diff
			Scat_limit=max_obs_diff[0]*(0.5)
 
 			if(scat_limit>2*(Scat_SD[i]/scat_avg[i]))
 				scat_limit=2*Scat_SD[i]/scat_avg[i]
 			endif
 
 
			Diff_calc_k_wave(possible_abs,k_range)
			wave max_obs_diff
			abs_limit=max_obs_diff[0]*(0.5)
 
  			if(abs_limit>2*(abs_SD[i]/abs_Avg[i]))
 				abs_limit=2*abs_SD[i]/abs_Avg[i]
 			endif
			 
			
			
			
			
			
			For (h=0; h<numpts_RI;h+=1) //cycles though each RI input point
						
			 
			
				if (CombinedBarrow[h][i]<RI_Calc(i))  //checks to see if the difference is the lowest one for that time point
						if(WheelBarrow_abs[h][i]/abs_avg[i]<abs_limit)
							if(WheelBarrow[h][i]/scat_avg[i]<scat_limit)
						
							RI_Calc[i]=combinedbarrow[h][i]  // if the difference is the lowest yet found, the dummy wave is set to that value
							RI_Results[i]=RI_Range(h)  // if the difference is the lowest yet found for timepoint , the n results is set to the corresponding n input value
							K_Results[i]=K_Range(h)  // if the difference is the lowest yet found timepoint, the k results is set to the corresponding k input value
							endif
						endif
				endif
			endfor
		endfor	

	
			duplicate/O RI_results,results_n
			duplicate/o K_Results, results_k

	killwaves RI_results, K_results, abs_Series, Back, combinedbarrow, conversion, disposable, disposable_abs, dn, ext, output, output_scat, output_abs, Qe_,RI_Calc,Scat_avg, Scat_Series,template,trowel,trowel_abs,wheelbarrow,wheelbarrow_abs,absolutebarrow,absolutebarrow_abs, k_results

	
	End




////simulates the theoretical extinction t series for a given time series of refractive index values and size distributions.	
	
	
function Simulate_ext_t_series(n_t_series,k_t_series,size_dist,wavelength, dia_bins)
wave n_T_series, k_t_series, size_dist,dia_bins
variable wavelength


variable i=0,p=0

make /o /n=1 n_range, k_range
duplicate /o n_t_series, scat_coef
duplicate /o k_t_series, abs_coef

for (i=0;i<numpnts(n_t_series);i+=1)  //runs through each RI point
	
n_range=n_t_series[i]
k_range=k_t_series[i]


efficency_matrix(n_range, k_range, wavelength, Dia_Bins) //calculates the efficiencies at every bin diameter for that given time series point

wave Eff_Matrix_k,Eff_Matrix_n

duplicate /o dia_bins, size_profile
for (p=0;p<numpnts(dia_bins);p+=1)
	
		size_profile[p]=size_dist[p][i]
endfor


	Coef_Calc(dia_bins,Eff_Matrix_n, Eff_matrix_k,Size_profile)
	
wave ext_func_n
wave ext_func_k


scat_coef[i]=ext_func_n[0]
abs_coef[i]=ext_func_k[0]

endfor
end   
