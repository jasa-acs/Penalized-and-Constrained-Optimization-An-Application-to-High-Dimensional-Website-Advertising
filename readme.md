# Penalized and Constrained Optimization: An Application to High-Dimensional Website Advertising

# Author Contributions Checklist Form

## Data

### Abstract

The case study considered in this paper is implemented using a subset of the 2011 comScore Media Metrix data. comScore's data is a commercial, proprietary data set purchased through comScore and accessed through the Wharton Research Data Service (www.wrds.upenn.edu). comScore records daily webpage usage information from a panel of 100,000 Internet users, whose behavior is recorded anonymously by individual computer. Using these comScore by-computer records, we construct a matrix of all websites visited and the number of times each computer visited each website (and how many webpages were viewed at each visit) during a particular time period.

The data used in the case study utilizes website visits by the 100,000 comScore users during January 2011. To create a manageable dataset, we manually identify the 500 most-visited websites in January 2011 which also supported Internet display ads. Thus, our filtered data contains a record of every computer which visited at least one of the 500 most-visited websites at least once (48,628 users) during the month of January. The NCL case study ultimately uses a matrix of 48,628 comScore users by 500 websites, where the matrix entries are the total number of webpages viewed by each user at each website during the month of January.

### Availability
Due to the contract signed with comScore at the time of purchase, this data *unfortunately cannot be made publicly available*. This data set is proprietary, purchased through comScore, Inc. The contract heavily restricts even characteristics of the data (for example, information on specific websites that appear in the data set). If desired, the authors can provide a copy of the contract for verification that the data set is unable to be distributed.

Because of this, we have also provided numerous simulations in the paper using data we have generated ourselves. For verification purposes, we have provided the complete code to run these simulations and generate their results in their entirety, as this has been where the most interest has been generated from our paper thus far. In addition, we provide a minimal working example to demonstrate the type of results shown in the proprietary case study.

### Description

(See above for proprietary data description. The data used for the simulation results is generated through our given “generate.data” function and thus is generated through the code itself and has no restrictions.)


## Code

### Abstract

The code includes all functions necessary to run the results found in Section 5 of the paper (the results run on simulated, non-proprietary data). We encourage readers to see the README file accompanying the submitted data and code for examples of the reproducible results code, as well as full details on the functions (including arguments and outputs).

In addition, we have included a collection of code and data files, also summarized in the README file:

*	PaC_Functions_Revised.R
*	PaC_Full_Reproducibility.R
    * PaC_Simple_Table2.R
    * PaC_Tables_2_and_3.R
    * PaC_Timings_Fig2.R
*	Minimal_Working_Example_Revised.R
*	Page_View_Matrix_Example.csv
*	500_Site_Info_Example.csv

The complete file for running the results in Section 5 is PaC_Full_Reproducibility.R. To get the exact results in the paper, you should be able to run this entire file directly. It sources the main function code, PaC_Functions_Tested_Revision.R, so the function R file should be stored in the main directory in order to run the code successfully.

However, because this complete file is lengthy and difficult to parse, we have included three subfiles for particular parts of the reproducibility depending on which part you wish to reproduce. Note that the code for Tables 2 and 3 is extremely similar, so we include PaC_Simple_Table2.R for users who wish to see a fully-commented version of the type of code used to generate the Tables. It will reproduce only the first line of Table 2, but should be much easier to follow through for understanding and usability. In addition, it should take approximately one minute to run complete for the full 100 iterations used in the paper.

To supplement this, we also include PaC_Tables_2_and_3.R, which is a file that can be run completely to reproduce the entirety of the Table 2 and Table 3 results. This file takes significantly longer to run (several hours for the full 100 iterations), but it will demonstrate the complete results of both Tables.

Finally, we also include PaC_Timings_Fig2.R. This file contains the code necessary to recreate a version of Figure 2 in Section 5.3 of the PaC paper. Because this code takes considerable time to run due to the changing values of the number of coefficients, this code is separated out for easier walkthrough. Please note that your run of Figure 2 will likely look similar but not exact to the presented Figure 2 in the PaC paper. This is because the complete timings code to generate Figure 2 is run several times and averaged, along with increasing maximum iterations and changing step size as necessary. Since this process takes several days to complete, the code provided here is done for a single run with a low default maximum number of iterations (12) and relatively large default step size (0.2). Even so, this code will take several hours to run. You may wish to modify the Figure 2 code to get an even more reasonable runtime, though this will further separate the generated results from Figure 2.

Last but not least, we have provided a minimal working example to demonstrate to users how to use the PaC method in the case study example. Because the data set in the case study is proprietary, we have also provided a CSV (Page_View_Matrix_Example.csv) file that shows an example of the structure of the data used in the paper. It is a 5000-user by 500-website matrix that contains data generated to mimic the comScore page view data. It is used in conjunction with an example of the information gathered regarding individual websites (500_Site_Info_Example.csv). These two files can be used with the Minimal_Working_Example_Revised.R code file to create examples of the results, as shown at the end of this file. (See Reproducibility in Instructions for Use below.)

The PaC_Functions_Revision.R file contains additional minimal functions for working with this data, as well as a binomial implementation of the main PAC function as well (as the case study implementation is based on a binomial function).

### Description

This code is delivered via the files described above. There is no particular licensing information for this code (and it will be developed into a complete package after the submission of this paper). The code will be hosted on authors’ websites as well as through the CRAN repository once the package is completed and the paper is accepted.

### Optional Information

No general hardware is required for this code.

Software requirements:

R packages to run PAC functions:

*	lars library (suggested version 1.2)
*	limSolve library (suggested version 1.5.5.3)
*	quadprog library (suggested version 1.5-5)

R packages to run reproducible code:

*	MASS library (7.3-32 or later)
*	MBESS library (4.0.0 or later)

In addition, the authors recommend the use of the 64-bit version of R (version 3.1.0 or later) if running the files in their entirety due to memory allocation, though this is not necessary.


## Instructions for use

### Reproducibility 

Tables 2 and 3 can be reproduced in their entirety, as well as Figure 2. In addition, we have supplied a sample of data generated to mimic the proprietary data (Page_View_Matrix_Example.csv) as well as example cost information by website (500_Site_Info_Example.csv). These can be used with the Minimal_Working_Example.R file to reproduce the example.

The main reproducible contribution without the proprietary data is the more general simulation results of Section 5, which could be implemented on any data set. We describe the workflow to reproduce these below, but encourage readers to refer to the README file for a more complete description.

To reproduce the analyses of Section 5 and the plot above, the file PaC_Full_Reproducibility.R can be run directly. Again, however, the authors recommend you run the subfile for whichever results you wish to reproduce rather than running the full file for reproducibility. Whichever file you wish to run, however, you will need to install the following packages if not already installed:

*	MBESS (4.0.0 or later) 
*	MASS (7.3-32 or later)
*	limSolve (suggested version 1.5.5.3)
*	lars (suggested version 1.2)
*	quadprog (suggested version 1.5-5)

The following workflow should be followed if running the full PaC_Full_Reproducibility.R code, if not running the code in its entirety at once:

1.	Install and load necessary packages
2.	Source the function code, PaC_Functions_Revised.R
3.	Set number of iterations using iter (default is 100, which will exactly reproduce the results of the paper but may take considerable time to run): line 27
4.	Run the sections of the Reproducibility Code file in order. They are:
    a.	Coefficient Path Setup (running the PAC fits for Table 3): line 33
    b.	Coefficient Path Setup with Error (running the PAC fits for Table 3): line 155
    c.	Lasso and Relaxed Fit Paths, used in Table 2: line 247
    d.	Calculation of the Error Rates in Table 2: line 1169
    e.	Lasso and Relaxed Fit Paths, used in Table 3: line 1238
    f.	Calculation of the Error Rates in Table 3: line 3027
    g.	Comparing Calculation Times (Figure 2): line 3176

## Notes

Because the data is proprietary and not publicly available, it is necessary to purchase the data through the comScore company directly. This process usually begins online through their Analytics portal (https://www.comscore.com/Products/Audience-Analytics/Media-Metrix), but you can also contact a representative directly through the Wharton Research Data Services (https://wrds-web.wharton.upenn.edu/wrds/index.cfm). The data is generally made available via queries on WRDS once access is granted, though comScore was able to provide CSV samples of the data as well upon request.
