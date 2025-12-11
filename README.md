# This repository contains Matlab codes used to analyze X-ray diffraction (XRD) synchrotron data collected at Soleil during the measurement week between 6 July 2021 to 12 July 2021.

To use these 4 codes, place the desired file(s) in a folder together with the corresponding dataset and the file "backcor.m".

Datasets can be found here: LINK

Modify the "scan_temp" variable at line 6 to analyze another dataset. This value must be equal to the one in the dataset file name.

1 - "Code_Pt_Data_Treatment.m" : This code performs peak fitting analysis on XRD spectra obtained from the Pt sample at a specified temperature. For each diffractogram, the code fits the peaks and calculates the Pt lattice parameter and the sample displacement. 

Available datasets for this code:
- "Data Pt - 20C.mat" (use scan_temp = '20')
- "Data Pt - 50C.mat" (use scan_temp = '50')
- "Data Pt - 100C.mat" (use scan_temp = '100')
- "Data Pt - 200C.mat" (use scan_temp = '200')
- "Data Pt - 250C.mat" (use scan_temp = '250')
- "Data Pt - 300C.mat" (use scan_temp = '300')
- "Data Pt - 350C.mat" (use scan_temp = '350')
- "Data Pt - 400C.mat" (use scan_temp = '400')
- "Data Pt - 450C.mat" (use scan_temp = '450')
- "Data Pt - 500C.mat" (use scan_temp = '500')
- "Data Pt - 550C.mat" (use scan_temp = '550')
- "Data Pt - 600C.mat" (use scan_temp = '600')
- "Data Pt - 650C.mat" (use scan_temp = '650')
- "Data Pt - 700C.mat" (use scan_temp = '700')
- "Data Pt - 735C.mat" (use scan_temp = '730')

2 - "Code_LSCFM_RT_Data_Treatment.m" : This code performs peak fitting analysis on XRD spectra obtained from the combinatorial LSCFM sample at room temperature (RT). For each diffractogram, the code fits the peaks and calculates the LSCFM lattice parameter using the mean sample displacement values from Pt measurements. Therefore, this code also requires the presence of the "zero_shift_mean.mat" file in the same folder as the code.

Available datasets for this code:
- "Data LSCFM - 20C.mat" (use scan_temp = '20')

3 - "Code_LSCFM_HT_Data_Treatment.m" : This code performs peak fitting analysis on XRD spectra obtained from the combinatorial LSCFM sample at a specified high temperature (HT). For each diffractogram, the code fits the peaks and calculates the LSCFM lattice parameter using the mean sample displacement values from Pt measurements. Therefore, this code also requires the presence of the "zero_shift_mean.mat" file in the same folder as the code.

Available datasets for this code:
- "Data LSCFM - 100C.mat" (use scan_temp = '100')
- "Data LSCFM - 200C.mat" (use scan_temp = '200')
- "Data LSCFM - 300C.mat" (use scan_temp = '300')
- "Data LSCFM - 400C.mat" (use scan_temp = '400')
- "Data LSCFM - 500C.mat" (use scan_temp = '500')
- "Data LSCFM - 600C.mat" (use scan_temp = '600')
- "Data LSCFM - 700C.mat" (use scan_temp = '700')
- "Data LSCFM - 735C.mat" (use scan_temp = '735')

4 - "Code_Pt_LSCFM_Data_Treatment.m" : This code performs peak fitting analysis on XRD spectra obtained from the combinatorial LSCFM sample with Pt at a specified high temperature (HT). For each diffractogram, the code fits the peaks and calculates the Pt lattice parameter and the sample displacement, or the LSCFM lattice parameter using the mean sample displacement values from Pt measurements. Therefore, this code also requires the presence of the "zero_shift_mean.mat" file in the same folder as the code.

Available datasets for this code:
- "Data Pt-LSCFM - 400C.mat" (use scan_temp = '400')
- "Data Pt-LSCFM - 700C.mat" (use scan_temp = '700')
