~ Installing RMySQL on Windows:
    1) Install latest RTools/cygwin
    2) Add R bin directory to path (right-click My Computer >> properties >> advanced system settings >> Environment Variables >> edit PATH append ";C:\Program Files\R\R-2.14.0\bin\x64\" sans quotes)
    3) Install MySQL-server
    4) Create or edit file C:\Program Files\R\R-2.14.0\etc\Renviron.site and add line like MYSQL_HOME=C:/Program Files/MySQL/MySQL Server 5.5 (path to your mysql files)
    5) Copy libmysql.lib from C:\Program Files\MySQL\MySQL Server 5.5\lib to C:\Program Files\MySQL\MySQL Server 5.5\lib\opt (silly, but thats where it looks...to meet dependencies.)
    6) Copy libmysql.dll to C:\Program Files\R\R-2.14.0\bin
    7) Copy libmysql.dll to C:\Program Files\MySQL\MySQL Server 5.5\bin\ (again, silly, but thats where it looks...)
    8) Run install.packages('RMySQL',type='source') and wait while compilation will end.

~ Running the parallel version of getNearestFeature, getSitesInFeature, getFeatureCounts:
    1) Load one of the following libraries depending on machine/OS: doParallel, doMC, doSMP, doSNOW, doMPI. 
    
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
        
        Example 4: library(doParallel); 
                   cl <- makeCluster(2); 
                   registerDoParallel(cl); 
                   getNearestFeature(..., parallel=T)  
                   
    3) Few backends launch workers in the background, so be sure to close them. Read the documentation of respective do* package to get more information. Few examples are shown below.
        For doSMP library, use stopWorkers(w)
        For doSNOW & doParallel library, use stopCluster(cl)
        