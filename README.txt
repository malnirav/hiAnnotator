~ Running the parallel version of getNearestFeature, getSitesInFeature, getFeatureCounts:
    1) Load one of the following libraries depending on machine/OS: doMC, doSMP, doSNOW, doMPI. 
    
    2) Register the parallel backend using registerDoXXXX() function depending on the library. See the examples below:  
        Example 1: library(doSMP); 
                   w <- startWorkers(2); 
                   registerDoSMP(w); 
                   getNearestFeature(..., parallel=T)        
                   
        Example 2: library(doMC); 
                   registerDoMC(2); 
                   getNearestFeature(..., parallel=T)        
                   
        Example 3: library(doSNOW); 
                   cl <- makeCluster(2, type = "SOCK"); 
                   registerDoSNOW(cl); 
                   getNearestFeature(..., parallel=T)                   
                   
    3) Few backends launch workers in the background, so be sure to close them. Read the documentation of respective do* package to get more information. Few examples are shown below.
        For doSMP library, use stopWorkers(w)
        For doSNOW library, use stopCluster(cl)
        