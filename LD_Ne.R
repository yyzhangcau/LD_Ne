beta_estimate<-function(LD_dataframe,alpha=1){
    df=LD_dataframe[LD_dataframe[,"r2"]>=0,]
    dist<-df[,"dist"]
    r2<-df[,"r2"]
    ls_model<-nls(r2~1/(alpha+beta*dist),start=list(beta=0),control = list(maxiter = 5000))
    return(coef(ls_model)['beta'])
}

calculate_c_for_each_window<-function(LD_matrix,window_size=500000,alpha=1,linkage_map=linkage_map,chr=1){
    #chromosome Length
    physical_L=max(LD_matrix[,"L2_pos"])
    #overlap=1/2 * window_size
    slides1=seq(0,physical_L,0.5*window_size)
    slides2=slides1+window_size
    #number of slides
    n_slides=length(slides1)
    names_rows=paste0(slides1/1000000,"-",slides2/1000000)
    out=data.frame(matrix(0,ncol=4,nrow=n_slides))
    colnames(out)<-c("physical_distances","N_SNPpairs","beta","c_estimates")
    out[,1]=names_rows
    ##calculate beta of the entire chromosome
    global_beta=beta_estimate(LD_matrix,alpha)
    ##calculate local beta
    for (i in 1:n_slides){
        index<-(LD_matrix[,"L1_pos"]>slides1[i])&(LD_matrix[,"L2_pos"]<slides2[i])
        df<-LD_matrix[index,]
        n_snppairs=dim(df)[1]
        out[i,2]=n_snppairs
        if(n_snppairs>=10){ # if number of SNP pairs less than 20, estimates may be not accurate. 
            beta=beta_estimate(df,alpha)
            out[i,3]=beta
        }
        else{
            out[i,3]=global_beta
            }
    }
    x=sum(out[,2]*out[,3]*window_size/1000000)*n_slides/(linkage_map[chr,2]*sum(out[,2]))
    out[,4]=out[,3]/x
    return(out)
    #has columns[1]physical_distances,[2]N_SNPpairs,[3]beta(coordinate with cM),[4]c_estimates(cM)
}

recombination_rate_accumulated_map<-function(df=out_c,window_size=500000){
    n_window=length(df[,1])+2
    physical_dist=seq(from=0,length.out=n_window,by=window_size/2)
    accumulated_c<-0

    for(i in 2:n_window-1){
        accumulated_c[i]<-sum(df[1:i-1,4])/2
    }
    accumulated_c[n_window]=accumulated_c[n_window-1]+df[n_window-1,4]/2
    c_map<-cbind(physical_dist,accumulated_c)
    return(c_map)
}

cM_for_each_snp<-function(pos,accumulated_map=c_map,window_size=500000){
    index=pos%/%(window_size/2)
    res=pos%%(window_size/2)
    pos_c<-accumulated_map[index+1,2]+(accumulated_map[index+2,2]-accumulated_map[index+1,2])/(window_size/2)*res
    return(pos_c)
}

calculate_c_for_each_pair<-function(LD_df,accumulated_map=c_map,alpha=1,window_size=500000){
    LD_df[,"L1_cM"]<-cM_for_each_snp(LD_df[,'L1_pos'],accumulated_map,window_size)
    LD_df[,"L2_cM"]<-cM_for_each_snp(LD_df[,'L2_pos'],accumulated_map,window_size)
    LD_df[,"dist"]=(LD_df[,'L2_cM']-LD_df[,"L1_cM"])/100
    return(LD_df)
}