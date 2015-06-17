
umi.data.location   <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\clean_data\\final_result_14-05-15-1551_2_mismatches.csv'
umi.data = read.csv(umi.data.location,row.names = 1)
umi.data.pnas <- pnas.umi.correction(umi.data,total.umis)
pnas.correction.location <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\clean_data\\pnas_correction_14-05-15-1551_2_mismatches.csv'

write.csv(umi.data.pnas,pnas.correction.location)




pnas.umi.correction <- function(k,m){
# m - the number of UMIs
# k  - the the number of unique labels that have been captured by  the design
    n <-   - m * log(1- k/m )  
    return(floor(n))
}
pnas.umi.correction <- Vectorize(pnas.umi.correction,'k')

total.umis <- 4^9

 write.csv
