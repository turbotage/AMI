Script order of H-scan on datasets 2 and 3

To compute scripts of H-scan on datasets 2 and 3 proceed according to the 
following instruction:

1. The main scripts for H-scan development are dataset2r.m and dataset3r.m. 
Run those first, but don't forget to change the filepath before.

2. At the end of the main scripts there is a section called "Additional 
analytical imagery" that is not included in the scope of our report but 
was utilized to better understand the signals spatial and temporal intensity 
propagation and was used to model the at the start- and the end included 
"Signal profile plots."

3. The final main script is called stdcomph.m (standard deviation comparation
 of h-scan), it computes the post-band pass filtering for isolating the 
contractions utlizing mybandpass.m as well as the performance measurement 
of contrast by calling the function roi_metric.m for dataset 2 before 
proceeding to its final task of using the H-scan adapted draw_std2.m script
to visualize the standard deviation mask of the chosen dataset. 
Note! You can run this script for either dataset 2 or 3 at a time. If you 
want to run them in succession you need to change the string chosen_dataset 
prior to the following run.

Addendum: sometimes I have varied climc values across some of the scripts so
or the exceptions in the functions written will not be caught, this is some-
thing I did not have time to tend to so just fill it in if it doesn't work.

/Rebecca Viklund