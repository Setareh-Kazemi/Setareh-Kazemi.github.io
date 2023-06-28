# This file contains the custom functions for the functional anova for fatigue assessment 


# * Read RPE Data ---------------------------------------------------------

read_rpe_data = function(file_name){
  
  # subject number consists of 4 lower case chrs (subj) and two digits
  subject_num = stringr::str_extract(file_name, pattern = "[:alnum:]{6}")
  
  # session number consists of 7 lower case chrs (session) and two digits
  session_num = stringr::str_extract(file_name, pattern = "[:alnum:]{9}")
  
  # returned rpe data
  readxl::read_excel(file_name) %>% 
    dplyr::transmute(
      subject = subject_num, 
      session = session_num,
      period = stringr::str_extract(type, '[:digit:]{1,}'), 
      time = time,
      rpe = RPE, 
      type = stringr::str_replace(type, '[:digit:]{1,}', '')
    ) %>% tibble::tibble()
  
}


# * Read Strength Data ----------------------------------------------------

read_strength_data = function(file_name){
  
  # subject number consists of 4 lower case chrs (subj) and two digits
  subject_num = stringr::str_extract(file_name, pattern = "[:alnum:]{6}")
  
  # session number consists of 7 lower case chrs (session) and two digits
  session_num = stringr::str_extract(file_name, pattern = "[:alnum:]{9}")
  
  # period is one of '1st_45min|2nd_45min|3rd_45min' and num is first digit
  period_num = stringr::str_extract(file_name, pattern = "1st_45min|2nd_45min|3rd_45min") %>% 
    str_extract(pattern = '[:digit:]{1}')
  
  # strength measurement number is the first digit in the last string after a /
  text_separated_by_slash = stringr::str_split(file_name, pattern = '/', simplify = T)
  strength_measurement_num = 
    stringr::str_extract(
      text_separated_by_slash[length(text_separated_by_slash)], 
      pattern = '[:digit:]{1}')
  
  # maximum absolute load
  max_abs_load = readxl::read_excel(path = file_name) %>% 
    dplyr::pull(var = 2)  %>% min(na.rm = T) %>% abs()
  
  # returned values
  tibble::tibble(subject = subject_num, session = session_num, period = period_num,
                 strength_meas_num = strength_measurement_num, max_load = max_abs_load)
}



# * Replace NA for Last RPE by its previous value -------------------------

impute_time_last_rpe = function(index = NULL, df = period1_df){
  df[index, 'time'] = df[(index-1), 'time']
  df[index, 'rpe'] = df[(index-1), 'rpe']
  return(df)
}



# * Plotting Function -----------------------------------------------------

rpe_str_plot = function (subject_id = NULL, df = NULL, cols_to_plot = NULL, fig_path ='figs_SK/'){
  
  # converting our data frame to a tall data frame to facilitate the plotting of data
  df_tall = df %>% 
    # filtering to our specific id
    dplyr::filter(subject == subject_id) %>% 
    # converting to a tall dataset
    tidyr::pivot_longer(cols = dplyr::all_of(cols_to_plot), names_to = "response") %>% 
    # converting time to numeric, ordering levels for experimental condition
    dplyr::mutate(
      time = as.numeric(time),
      exp_condition = factor(
        exp_condition,
        levels = c('1.5-15', '2.5-5', '2.5-10', '2.5-15'))) %>% 
    # grouping by experimental condition and response so we can impute values
    # for better visualization of the trends in strength and rpe data
    dplyr::group_by(exp_condition, response) %>% 
    dplyr::mutate(
      value_imputted = imputeTS::na_interpolation(value, option = 'linear'))
  
  # creating a named vector of cleaned response names to be used in labeling
  # the facets created by the facet_grid function
  response_names = tools::toTitleCase(
    stringr::str_replace_all(unique(df_tall$response), pattern = '_', replacement = ' ')
  )
  
  response_vec = c(unique(df_tall$response), 'Response' )
  names(response_vec) = c(response_names, 'Statistic of Interest')
  
  # creating the ggplot
  rpe_str_plot = df_tall %>%
    # cleaning up the names for the response and exp_condition to make the graph
    # easier to interpret
    dplyr::mutate(
      response = tools::toTitleCase(stringr::str_replace_all(response, '_', ' ')) %>% 
        stringr::str_replace_all(pattern = 'rpe|Rpe', replacement = 'RPE'),
      exp_condition = 
        paste(stringr::str_replace_all(exp_condition, '-', 'kg - '), 'bpm')) %>% 
    # the ggplot
    ggplot2::ggplot(aes(x = time, group = response)) +
    ggplot2::geom_line(aes(y = value_imputted)) +
    ggplot2::geom_point(aes(y = value), size = 1.5) +
    ggplot2::facet_grid(rows = vars(response), cols = vars(exp_condition), scales = 'free_y') +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n=3), limits = c(0, 50)) +
    ggplot2::scale_color_brewer(palette = 'Dark2') +
    ggplot2::labs(
      x = 'Time', y = 'Statistic of Interest',
      title = paste('Strength and RPE Data for Subject',
                    stringr::str_extract(subject_id, pattern = '[:digit:]{2}')),
      subtitle = str_wrap('Dots show true values, and lines involve interpolation for visualization purposes', width=70),
      caption = str_wrap(paste0('The ', unique(df_tall$age), ' yr-old ', 
                       unique(df_tall$sex), 
                       ' subject had a BMI of ', unique(round(df_tall$bmi, 1)), 
                       ' with waist and hip circumfurances (cm) of ', 
                       unique(round(df_tall$waist_circumf_cm, 1)), ' and ', 
                       unique(round(df_tall$hip_circumf_cm, 1)), ', respectively.') , width = 57)
    ) +
    
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::theme(legend.position = 'none') +
    ggplot2::theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )+
    ggplot2::theme(
       plot.subtitle = ggplot2::element_text(size = 16),  # Increase subtitle size
       plot.caption = ggplot2::element_text(size = 16, hjust=0),  # Increase caption size
       axis.title = ggplot2::element_text(size = 14),  # Increase axis label size
       axis.text = ggplot2::element_text(size = 14) # Increase axis text size
    )
  
  # saving the plot
  ggplot2::ggsave(filename = paste0(fig_path, 'rpe_str_', subject_id, '.png'),
                  plot = rpe_str_plot,
                  dpi = 600, width = 6.5, height = 4.5)
}


# * Pivot Wider for ANOVA -------------------------------------------------

pivot_wider_anova = function(df = period1_df, response_var_name = 'frac_max_load', time_sep = 9){
  df_wide = df %>% 
    dplyr::select( subject, exp_condition, time, as.name(response_var_name) ) %>% 
    # converting to a wide df, where each col is combination of subject x exp condition
    tidyr::pivot_wider(
      names_from = c(subject, exp_condition), 
      values_from =as.name(response_var_name), 
      values_fn = max
    ) %>% 
    # keep time values divisible by time_sep since we measured response every time_sep minutes
    dplyr::filter(time %in% seq(0, 45, time_sep)) %>% 
    # making the time column our row names
    tibble::column_to_rownames(var = 'time') %>% 
    # selecting only the columns that do not contain any NA values
    dplyr::select( tidyselect::where(~ !any(is.na(.x)) ) )
  
  return(df_wide)
}

# * Multivariate Box Cox  -------------------------------------------------

# The "fam" function, is actually the "bcPower" function from the "powertransfomr" library, 
# is defined to calculate the normalized transformed data ("out" as the output), given the box-cox parameter lambda and the vector of the data. 
# First, "z" is calculated as the box cox transformation.
# Next, the normalized transformation is calculated as the product of the transformed data "z" and the geometric mean.

fam <- function(U, lambda, jacobian.adjusted=FALSE, gamma=NULL) {
  if(!is.null(gamma)) bcPower(t(t(as.matrix(U) + gamma)), lambda, jacobian.adjusted) else{
    bc1 <- function(U, lambda){
      if(any(U[!is.na(U)] <= 0)) stop("First argument must be strictly positive.")
      z <- if (abs(lambda) <= 1.e-6) log(U) else ((U^lambda) - 1)/lambda
      if (jacobian.adjusted == TRUE) {
        z * (exp(mean(log(U), na.rm=TRUE)))^(1-lambda)} else z
    }
    out <- U
    out <- if(is.matrix(out) | is.data.frame(out)){
      if(is.null(colnames(out))) colnames(out) <-
          paste("Z", 1:dim(out)[2], sep="")
      for (j in 1:ncol(out)) {out[, j] <- bc1(out[, j], lambda[j]) }
      colnames(out) <- paste(colnames(out), round(lambda, 2), sep="^")
      out}  else
        bc1(out, lambda)
    out}}

# The "estimateTransform.default" function calculates the maximum likelihood estimate of the box cox transformation parameter (lambda).The only modification on this function is adding the lower and upper limit for the estimated lambda when "res" is calculated.
# The output of this function is calculated based on the third formula provided in this article: https://journals.sagepub.com/doi/epdf/10.1177/1536867X1001000108

estimateTransform.default <- function(X, Y, weights=NULL,
                                      family="bcPower", start=NULL, method="L-BFGS-B", ...) { 
  fam <- "bcPower" 
  Y <- as.matrix(Y) # coerces Y to be a matrix.
  X <- as.matrix(X)# coerces X to be a matrix.
  w <- if(is.null(weights)) 1 else sqrt(weights)
  nc <- dim(Y)[2]
  nr <- nrow(Y)
  xqr <- qr(w * X)
  llik <- function(lambda){
    (nr/2)*log(((nr - 1)/nr) *
                 det(var(qr.resid(xqr, w*fam(Y, lambda, j=TRUE)))))
  }
  llik1d <- function(lambda,Y){
    (nr/2)*log(((nr - 1)/nr) * var(qr.resid(xqr, w*fam(Y, lambda, j=TRUE))))
  }
  if (is.null(start)) {
    start <- rep(1, nc)
    for (j in 1:nc){
      res<- suppressWarnings(optimize(
        f = function(lambda) llik1d(lambda,Y[ , j, drop=FALSE]),
        lower=-3, upper=+3))
      start[j] <- res$minimum
    }
  }
  res <- optim(start, llik, hessian=TRUE, method="L-BFGS-B", lower = -3, upper =3)
  if(res$convergence != 0)
    warning(paste("Convergence failure: return code =", res$convergence))
  res$start<-start
  res$lambda <- res$par
  names(res$lambda) <-
    if (is.null(colnames(Y))) paste("Y", 1:dim(Y)[2], sep="")
  else colnames(Y)
  roundlam <- res$lambda
  stderr <- sqrt(diag(solve(res$hessian)))
  lamL <- roundlam - 1.96 * stderr
  lamU <- roundlam + 1.96 * stderr
  for (val in rev(c(1, 0, -1, .5, .33, -.5, -.33, 2, -2))) {
    sel <- lamL <= val & val <= lamU
    roundlam[sel] <- val
  }
  res$roundlam <- roundlam
  res$invHess <- solve(res$hessian)
  res$llik <- res$value
  res$par <- NULL
  res$family<-family
  res$xqr <- xqr
  res$y <- Y
  res$x <- as.matrix(X)
  res$weights <- weights
  res$family<-family
  class(res) <- "powerTransform"
  res
}



conc_lik = function(z){
  YR = t(z)
  XR = rep(1, dim(YR)[1])
  transformationz = estimateTransform.default(XR,YR)
  lambdaz = transformationz$lambda
  z_T_val = car::bcPower(YR,lambdaz)
  
  list = list(lambdas = lambdaz,
              transformed_variable = z_T_val)
  return(list)
}

# * FANOVA --------------------------------------------------------------

#The fanova.RPm.bootS function below is a modified version of the "fanova.RPm.boot" function from the "fda.usc" library.
# The only modification is changing "aov" function to "Anova" to perform a type-II anova test instead of type-I.
# This modification leads to having consistent outputs if the order of factor in the given formula changes, and it will show the same p-values as the t-test of coefficients in a regression analysis.
# Another benefit of using this "Anova" function is showing the interaction effects in the summary of the results even if their df is equal to zero.

fanova.RPm.bootS=function(object,formula,data.fac,RP=RP,
                          alpha=alpha,zproj=NULL,par.zproj=list(norm=TRUE),nboot=nboot,contrast=NULL,pr=FALSE,...){
  if (is.data.frame(object)) object=as.matrix(object)
  else if (is.fdata(object)) object=object[["data"]]
  min.data.fac <- min(table(data.fac))
  # if (min.data.fac==0)  warning("Contingency table of factor levels (data.fac argument) contains 0 counts  values")
  nrow=nrow(object);ncol=ncol(object)
  nprRP=max(RP)
  # if (!is.null(zproj)) nprRP=nrow(zproj)
  if (!is.null(zproj)) {
    if (is.fdata(zproj)) { 
      if (nrow(zproj)>=nprRP) { 
        z=zproj[1:nprRP] 
      } else {
        stop(paste("Not enough functions in zproj",length(zproj)," to compute ",nprRP," projections"))  
      }
    } else if (is.matrix(zproj)){
      if (nrow(zproj)>=nprRP){
        z=zproj[1:nprRP,] } else {
          stop(paste("Not enough rows in zproj",nrow(zproj)," to compute ",nprRP," projections"))
        }
    }  else if (is.function(zproj)){
      if (is.fdata(object)) z=do.call(zproj,modifyList(list(n=nprRP,t=object$argvals),par.zproj))
      if (is.matrix(object)) z=do.call(zproj,modifyList(list(n=nprRP),par.zproj))
    } else {stop("Parameter zproj is neither an fdata object or a function")}
    
    #        modulo=function(z){sqrt(sum(z^2))}
    #        z=rnorm(ncol*nprRP)
    #        z=matrix(z,nrow=nprRP,ncol=ncol)
    #        modu=apply(z,1,modulo)
    #        z=z/modu
  }
  terms.fd=attr(terms(formula),"term.labels")
  nterms=length(terms.fd)+1
  lterms=length(terms.fd)
  ff=attr(terms.fd,"factors")
  ffcol=colnames(ff)
  fml=as.formula(paste("object ~ ", paste(terms.fd, collapse= "+")))
  fmlb=as.formula(paste("value ~ ", paste(terms.fd, collapse= "+")))
  if (is.null(contrast))  {
    bb=array(NA,dim=c(nboot,nprRP,nterms-1))
    if (lterms==1)    dimnames(bb)[[3]]=list(terms.fd)
    else     dimnames(bb)[[3]]=(terms.fd) #
    ncontrast=0
  }
  else {
    b=length(contrast)
    mf2=model.frame(formula,data.fac)
    uniq=apply(mf2,2,function(x){length(unique(x))})
    ncontrast=rep(0,len=b)
    name.contr=rep(NA,len=sum(ncontrast))
    bb2=length(uniq)
    cnombres=ncontrast=rep(0,len=b)
    tnombres=rep(FALSE,len=length(ncontrast))
    tgroups=rep(FALSE,len=(length(ffcol)+1))
    for (i in 1:b)    {
      a=which(names(contrast)[i]==colnames(mf2))
      tgroups[a]=tnombres[a]=TRUE
      if (is.vector(contrast[[i]])) {
        ncontrast[a]=1
        contrast[[i]]=matrix(contrast[[i]],ncol=1)
      }
      else ncontrast[a]=ncol(contrast[[i]])
      names(ncontrast)[a]=names(contrast[i])
    }
    j=1;ji=1;jk=1
    for (i in 1:length(ncontrast))    {
      if (tgroups[ji])  {
        name.contr[j:(j+ncontrast[i]-1)]=paste("C",j:(j+ncontrast[i]-1),".",
                                               names(contrast[jk]),sep="")
        colnames(contrast[[jk]])=name.contr[j:(j+ncontrast[i]-1)]
        j=j+ncontrast[i];jk=jk+1
      }
      ji=ji+1
    }
    mat=matrix(NA,ncol=(nterms+sum(ncontrast)-1),nrow=nprRP)
    colnames(mat)=c(terms.fd,name.contr)
    if (pr) {
      print(ncontrast);print("ncontrast")
      print(contrast);print("contrast")
      print(name.contr);print("name.contr")
    }
    bb=array(NA,dim=c(nboot,nprRP,(nterms+sum(ncontrast)-1)))
  }
  for (k in 1:ncol(data.fac)) assign(names(data.fac)[k],as.factor(data.fac[,k]))
  fit=manova(fml)
  for (i in 1:nboot){
    
    l=sample(1:nrow(fit$residuals),replace=TRUE)
    funcboot=fit$residuals[l,]
    for (j in 1:nprRP){
      if (is.fdata(z)) value=funcboot%*%z$data[j,]
      if (is.matrix(z)) value=funcboot%*%z[j,]
      #       value=funcboot%*%z[j,]
      
      mdata=as.data.frame(cbind(value,data.fac))
      colnames(mdata)=c("value",colnames(data.fac))
      mod <- lm( fmlb, mdata )
      out = Anova( mod )
      #mat[j,1:(nterms-1)]=Out[1:(nterms-1),4]
      
      #result=aov(fmlb,data=mdata)
      #out=summary(result)
      if (pr) print(out)
      if (!is.null(contrast)) {
        result3=lm(fmlb,data=mdata,contrasts=contrast,...)
        out4=summary(result3)
        if (pr)     print(out4)
        if (length(out)==1) {
          bb[i,j,1:(nterms-1)]=out[[1]][1:(nterms-1),5] #p-valor
          if (pr) print(out4$coefficients)
          ind=nterms
          ind2=2
          for (ib in 1:bb2)    {
            if (ncontrast[ib]!=0) {
              bb[i,j,ind:(ind+ncontrast[ib]-1)]=out4$coefficients[ind2:(ind2+ncontrast[ib]-1),4]
            }
            ind=ind+ncontrast[ib];ind2=ind2+ncontrast[ib]-1
          }       }   }
      else   bb[i,j,1:(nterms-1)]=out[1:(nterms-1),4] #out[[1]][1:(nterms-1),5] #new
    }
  }
  resboot=vector("list",length(RP))
  for (k in 1:length(RP)){
    if (RP[k]==1)              resboot[[k]]=bb[,1,]
    else {
      if (lterms==1)       resboot[[k]]=apply(bb[,1:RP[k],],c(1),min)
      else   resboot[[k]]=apply(bb[,1:RP[k],],c(1,3),min)
    }
  }
  names(resboot)=paste("RP",RP,sep="")
  return(resboot)
}


# The fanova function below is a modified version of the "fanova.RPm" function from the "fda.usc" library. 
# The main modification is changing the "aov" function to "Anova" function for the same reasons as explained above.
# Another modification is the method of aggregating the bootstrap p-values, 
# which instead of using the "ecdf" we are calculating the mean of p-values as proposed in [].

fanova = function(nboot=nboot, RP, object,formula,data.fac,alpha=alpha,
                  zproj=NULL,par.zproj=list(norm=TRUE),
                  pr=FALSE,w=rep(1,ncol(object)),contrast=NULL){
  set.seed(2023)
  if (class(object)[1] %in% c("matrix","data.frame")) {dataM=as.matrix(object)}else if(is.fdata(object)) {dataM=object[["data"]]}
  lfac=unlist(lapply(data.fac,is.factor))
  min.data.fac<-min(table(data.fac[,lfac]))
  # if (min.data.fac==0)  warning("Contingency table of factor levels (data.fac argument) contains 0 counts  values")
  nrow=nrow(dataM);ncol=ncol(dataM)
  bonf=(1-alpha)/RP
  nprRP=max(RP)
  
  terms.fd=attr(terms(formula),"term.labels")
  
  fml=as.formula(paste("value ~ ", paste(terms.fd, collapse= "+")))
  if (is.null(nrow) || is.null(ncol)) stop("fdata must be a matrix")
  nterms=length(terms.fd)+1
  #
  projMV=function(n,m,w=rep(1,m),norm=TRUE){
    set.seed(2023)
    modulo=function(z){sqrt(sum(z^2))}
    z=rnorm(n*m)
    z=matrix(z,nrow=n,ncol=m)
    z=t(t(z)*w)
    if (norm) {
      modu=apply(z,1,modulo)
      z=z/modu}
  }
  if ((class(object) %in% c("matrix","data.frame")) & is.null(zproj)) zproj=projMV
  if ((is.fdata(object)) & is.null(zproj)) zproj=rproc2fdata
  if (is.fdata(zproj) & is.fdata(object)) {
    if (nrow(zproj)>=nprRP) {
      z=zproj[1:nprRP]
    } else {
      stop(paste("Not enough functions in zproj",length(zproj)," to compute ",nprRP," projections"))
    }
  } else if (is.matrix(zproj) & (class(object)[1] %in% c("matrix","data.frame"))){
    if (nrow(zproj)>=nprRP){
      z=zproj[1:nprRP,] } else {
        stop(paste("Not enough rows in zproj",nrow(zproj)," to compute ",nprRP," projections"))
      }
  }  else if (is.function(zproj)){
    if (is.fdata(object)) {z=do.call(zproj,modifyList(list(n=nprRP,t=object$argvals),par.zproj))}
    else if (class(object)[1] %in% c("matrix","data.frame")) {z=do.call(zproj,modifyList(list(n=nprRP,m=ncol(object)),par.zproj))}
  } else {stop("Parameter zproj is neither an fdata object or a function")}
  
  ff=attr(terms.fd,"factors")
  ffcol=colnames(ff)
  if (is.null(contrast)){
    mat=matrix(NA,ncol=(nterms-1),nrow=nprRP)
    colnames(mat)=c(terms.fd)
    ncontrast=0
  } else {
    b=length(contrast)
    mf2=model.frame(formula,data.fac)
    uniq=apply(mf2,2,function(x){length(unique(x))})
    ncontrast=rep(0,len=b)
    name.contr=rep(NA,len=sum(ncontrast))
    bb=length(uniq)
    cnombres=ncontrast=rep(0,len=b)
    tnombres=rep(FALSE,len=length(ncontrast))
    tgroups=rep(FALSE,len=(length(ffcol)+1))
    if (pr) {print(contrast);print("contrast     contrast  contrast")}
    for (i in 1:b)    {
      a=which(names(contrast)[i]==colnames(mf2))
      tgroups[a]=tnombres[a]=TRUE
      if (is.vector(contrast[[i]])) {
        if (pr) print("Is vector")
        print(a)
        ncontrast[a]=1
        contrast[[i]]=matrix(contrast[[i]],ncol=1)
      }
      else ncontrast[a]=ncol(contrast[[i]])
      names(ncontrast)[a]=names(contrast[i])
    }
    j=1;ji=1;jk=1
    for (i in 1:length(ncontrast))    {
      if (tgroups[ji])  {
        name.contr[j:(j+ncontrast[i]-1)]=paste("C",j:(j+ncontrast[i]-1),".",
                                               names(contrast[jk]),sep="")
        colnames(contrast[[jk]])=name.contr[j:(j+ncontrast[i]-1)]
        j=j+ncontrast[i];jk=jk+1
      }
      ji=ji+1
    }
    mat=matrix(NA,ncol=(nterms+sum(ncontrast)-1),nrow=nprRP)
    colnames(mat)=c(terms.fd,name.contr)
    if (pr) {
      print(ncontrast);print("ncontrast")
      print(contrast);print("contrast")
      print(name.contr);print("name.contr")
    }
  }
  
  library("car")
  
  for (j in 1:nprRP){
    #    value=data%*%z[j,]
    if (is(object,"fdata")) {
      value=inprod.fdata(object,z[j])
    } else if (is(object,"data.frame") | is(object,"matrix")) {
      value=dataM%*%z[j,]
    }
    mdata=as.data.frame(cbind(value,data.fac))
    colnames(mdata)=c("value",colnames(data.fac))
    
    mod <- lm( fml, mdata )
    Out = Anova( mod )
    mat[j,1:(nterms-1)]=Out[1:(nterms-1),4]
    
    #   result=aov(fml,data=mdata)
    #   out=summary(result)
    
  }
  if (pr) {print(summary(mat));print(paste("Bonferroni:",round(bonf,6)))}
  tmat=matrix(NA,nrow=length(RP),ncol=ncol(mat));colnames(tmat)=colnames(mat);rownames(tmat)=paste("RP",RP,sep="")
  mins=matrix(NA,nrow=length(RP),ncol=ncol(mat));colnames(mins)=colnames(mat);rownames(mins)=rownames(tmat)
  tFDR=matrix(NA,nrow=length(RP),ncol=ncol(mat));colnames(tFDR)=colnames(mat);rownames(tFDR)=rownames(tmat)
  pFDR=matrix(NA,nrow=length(RP),ncol=ncol(mat));colnames(pFDR)=colnames(mat);rownames(pFDR)=rownames(tmat)
  for (l in 1:length(RP)){
    if (RP[l]==1)      tmat[l,]=mat[1,]
    else {
      tmat[l,1:(nterms-1+sum(ncontrast))]=
        apply(matrix(mat[1:RP[l],1:(nterms-1+sum(ncontrast))],ncol=(nterms-1+sum(ncontrast))),2,min)}
    if (RP[l]==1) mins[l,]=rep(1,len=ncol(mat))
    else {
      try(mins[l,] <- apply(mat[1:RP[l],1:(nterms-1+sum(ncontrast))],2,which.min), silent = T)
    }
    if (RP[l]==1) tFDR[l,]=(mat[1,]<(1-alpha))
    else {
      tFDR[l,1:(nterms-1+sum(ncontrast))]=
        apply(matrix(mat[1:RP[l],1:(nterms-1+sum(ncontrast))],ncol=(nterms-1+sum(ncontrast))),2,FDR,alpha=alpha)}
    if (RP[l]==1) pFDR[l,]=mat[1,]
    else {
      pFDR[l,1:(nterms-1+sum(ncontrast))]=
        apply(matrix(mat[1:RP[l],1:(nterms-1+sum(ncontrast))],ncol=(nterms-1+sum(ncontrast))),2,pvalue.FDR)}
  }
  tbonf = (tmat < matrix(bonf, nrow = length(RP), ncol = ncol(mat)))
  colnames(tbonf) = colnames(mat);rownames(tbonf) = rownames(tmat)
  pbonf=tmat*RP
  colnames(pbonf) = colnames(mat);rownames(pbonf) = rownames(tmat)
  if (is.null(contrast)) {
    resb=matrix(NA,nrow=length(RP),ncol=nterms-1)    } else { resb=matrix(NA,nrow=length(RP),ncol=(nterms-1+sum(ncontrast)))}
  if (pr) {print(tbonf);print(pbonf);print(tFDR);print(pFDR)}
  if (nboot>0){
    if (pr) print("bootstrap procedure") ############
    resboot=fanova.RPm.bootS(dataM,formula,data.fac,RP=RP,alpha=alpha,
                             nboot=nboot,zproj=z,
                             contrast=contrast,pr=pr)
    if (pr) print(resboot)
    if (is.null(contrast)){nbuc=nterms-1}
    else {nbuc=nterms-1+sum(ncontrast)}
    for (l in 1:length(RP)){
      for (k in 1:nbuc){
        #if ((nterms-1)==1)    #Fn=ecdf(resboot[[l]])
        #else if(is.na(resboot[[l]][,k])) resb[l,k]<-NA
        #else #Fn=ecdf(resboot[[l]][,k])
        resb[l,k]=mean(resboot[[l]][,k])
        #Fn(tmat[l,k])
      }}
    rownames(resb)=paste("RP",RP,sep="")
  }
  pFDR[pFDR>1]=1
  pbonf[pbonf>1]=1
  res=list("proj"=z,"mins"=mins,"result"=mat,
           "test.Bonf"=tbonf,"p.Bonf"=pbonf,"test.FDR"=tFDR,"p.FDR"=pFDR, "resboot"=resboot, "p.Boot" = resb)
  if (nboot>0) {
    resb[resb>1]=1
    res$test.Boot=(resb<(1-alpha));res$p.Boot=resb
    colnames(res$p.Boot)=colnames(mat);colnames(res$test.Boot)=colnames(mat)}
  if (!is.null(contrast)) res$contrast <- contrast
  class(res) <- "fanova.RPm"
  return(res)
}


# * Smooth Function ---------------------------------------------------------

smooth_functional_var = function(nbasis, var){
  basis_var <- fda::create.bspline.basis(c(0,45), nbasis= nbasis)
  smooth_var <- fda::smooth.basis(argvals = seq(0,45, length.out = dim(t(var))[1]), y = t(var), fdParobj = basis_var)$fd
  list = list(basis = basis_var,
              smooth_var = smooth_var)
  return(list)
  }

# * Functional Clustering ---------------------------------------------------------

functional_clustering = function(var, k){
  functional_cluster = funHDDC::funHDDC(var,k,model=c('ABQkDk'), init='kmeans', threshold=0.2)
  list = list(functional_cluster = functional_cluster,
              cluster_label = functional_cluster$class)
  
}

# * Plot the Smooth Variables --------------------------------------------

plot_smooth = function(smooth_RPE, smooth_TRPE, smooth_STR, smooth_TSTR, cluster , smooth_factor = 0.2, fig_path ='figs_SK/'){
  # Changing fdobjects to dataframes to plot them with ggplot
  
  x_plot <- seq(0,45,smooth_factor)
  df_RPE <- data.frame(x=x_plot, yhat = predict(smooth_RPE, newdata=x_plot))
  df_TRPE <- data.frame(x=x_plot, yhat = predict(smooth_TRPE, newdata=x_plot))
  df_STR <- data.frame(x=x_plot, yhat = predict(smooth_STR, newdata=x_plot))
  df_TSTR <- data.frame(x=x_plot, yhat = predict(smooth_TSTR, newdata=x_plot))
  
  #creating long dataframes and adding the functional and gender cluster labels
  df_long_RPE <- df_RPE %>%
    pivot_longer(cols = starts_with("y"), names_to = "groups", values_to = "RPE") %>% 
    dplyr::mutate(c2 = factor(rep(functional_cluster_2_TRPE,length(x_plot))), 
                  cg2 = factor(rep(functional_cluster_g2_TRPE,length(x_plot))),
                  cg = factor(rep(period1_rpe_transformed_w_labels$sex,length(x_plot))))
  
  df_long_TRPE <- df_TRPE %>%
    pivot_longer(cols = starts_with("y"), names_to = "groups", values_to = "TRPE")%>% 
    dplyr::mutate(c2 = factor(rep(functional_cluster_2_TRPE,length(x_plot))), 
                  cg2 = factor(rep(functional_cluster_g2_TRPE,length(x_plot))),
                  cg = factor(rep(period1_rpe_transformed_w_labels$sex,length(x_plot))))
  
  df_long_STR <- df_STR %>%
    pivot_longer(cols = starts_with("y"), names_to = "groups", values_to = "STR")%>% 
    dplyr::mutate(c2 = factor(rep(functional_cluster_2_TSTR,length(x_plot))), 
                  cg2 = factor(rep(functional_cluster_g2_TSTR,length(x_plot))),
                  cg = factor(rep(period1_strength_transformed_w_labels$sex,length(x_plot))))
  
  df_long_TSTR <- df_TSTR %>%
    pivot_longer(cols = starts_with("y"), names_to = "groups", values_to = "TSTR")%>% 
    dplyr::mutate(c2 = factor(rep(functional_cluster_2_TSTR,length(x_plot))), 
                  cg2 = factor(rep(functional_cluster_g2_TSTR,length(x_plot))),
                  cg = factor(rep(period1_strength_transformed_w_labels$sex,length(x_plot))))
  
  #plots with clusters
  
  cluster <- sym(cluster)
  
  plot_RPE = ggplot(data=df_long_RPE, aes(x=x, y=RPE, group = groups, color = !!cluster)) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=5), limits = c(0,46))+
    scale_y_continuous(breaks=scales::pretty_breaks(n=5), limits = c(0,10))+
    labs(color = "Cluster")+
    geom_line(aes(x=x, y=RPE)) +
    theme_bw(base_size = 12)+
    scale_color_grey()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x="Time", y="RPE") 
  
  plot_TRPE = ggplot(data=df_long_TRPE, aes(x=x, y=TRPE, group = groups, color = !!cluster)) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=5), limits = c(0,46))+
    scale_y_continuous(breaks=scales::pretty_breaks(n=5), limits = c(0,10))+
    labs(color = "Cluster")+
    geom_line(aes(x=x, y=TRPE)) +
    theme_bw(base_size = 12)+
    scale_color_grey()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x="Time", y="TRPE") 
  
  plot_STR = ggplot(data=df_long_STR, aes(x=x, y=STR, group = groups, color = !!cluster)) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=5), limits = c(0,46))+
    scale_y_continuous(breaks=scales::pretty_breaks(n=5), limits = c(-0.5,1.5))+
    labs(color = "Cluster")+
    geom_line(aes(x=x, y=STR)) +
    theme_bw(base_size = 12)+
    scale_color_grey()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x="Time", y="STR")
  
  plot_TSTR = ggplot(data=df_long_TSTR, aes(x=x, y=TSTR, group = groups, color = !!cluster)) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=5), limits = c(0,46))+
    scale_y_continuous(breaks=scales::pretty_breaks(n=5), limits = c(-0.5,1.5))+
    labs(color = "Cluster")+
    geom_line(aes(x=x, y=TSTR)) +
    theme_bw(base_size = 12)+
    scale_color_grey()+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x="Time", y="TSTR")
  
  smooth_plot = ggpubr::ggarrange(
    plot_RPE , plot_TRPE, plot_STR, plot_TSTR,ncol = 2, nrow =2, common.legend = T
  )
  
  ggplot2::ggsave(filename = paste0(fig_path, 'plot_smooth_', rlang::as_string(cluster) , '.png'),
                  plot = smooth_plot,
                  dpi = 600, width = 6.5, height = 4.5)
  return(smooth_plot)
}

# * Plotting Function for One Subject -----------------------------------------------------

rpe_str_plot_one = function (subject_id = NULL, df = NULL, cols_to_plot = NULL, fig_path ='figs_SK/'){
  
  # converting our data frame to a tall data frame to facilitate the plotting of data
  df_tall = df %>% 
    # filtering to our specific id
    dplyr::filter(subject == subject_id) %>% 
    # converting to a tall dataset
    tidyr::pivot_longer(cols = dplyr::all_of(cols_to_plot), names_to = "response") %>% 
    # converting time to numeric, ordering levels for experimental condition
    dplyr::mutate(
      time = as.numeric(time),
      exp_condition = factor(
        exp_condition,
        levels = c('1.5-15', '2.5-5', '2.5-10', '2.5-15'))) %>% 
    # grouping by experimental condition and response so we can impute values
    # for better visualization of the trends in strength and rpe data
    dplyr::group_by(exp_condition, response) %>% 
    dplyr::mutate(
      value_imputted = imputeTS::na_interpolation(value, option = 'linear'))
  
  # creating a named vector of cleaned response names to be used in labeling
  # the facets created by the facet_grid function
  response_names = tools::toTitleCase(
    stringr::str_replace_all(unique(df_tall$response), pattern = '_', replacement = ' ')
  )
  
  response_vec = c(unique(df_tall$response), 'Response' )
  names(response_vec) = c(response_names, 'Statistic of Interest')
  
  # creating the ggplot
  rpe_str_plot_one = df_tall %>%
    # cleaning up the names for the response and exp_condition to make the graph
    # easier to interpret
    dplyr::mutate(
      response = tools::toTitleCase(stringr::str_replace_all(response, '_', ' ')) %>% 
        stringr::str_replace_all(pattern = 'rpe|Rpe', replacement = 'RPE'),
      exp_condition = 
        paste(stringr::str_replace_all(exp_condition, '-', 'kg - '), 'bpm')) %>% 
    # the ggplot
    ggplot2::ggplot(aes(x = time, group = response)) +
    ggplot2::geom_line(aes(y = value_imputted)) +
    ggplot2::geom_point(aes(y = value), size = 1.5) +
    ggplot2::facet_grid(rows = vars(response), cols = vars(exp_condition), scales = 'free_y') +
    # ggplot2::scale_y_continuous(
    #   breaks = scales::pretty_breaks(n = 5),
    #   limits = (ifelse(df_tall$response == "RPE", c(0, 10), c(0, 1)))
    # ) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n=3), limits = c(0, 50)) +
    ggplot2::scale_color_brewer(palette = 'Dark2') +
    ggplot2::labs(
      x = 'Time', y = 'Statistic of Interest',
      title = paste('Strength and RPE Data for Subject',
                    stringr::str_extract(subject_id, pattern = '[:digit:]{2}')),
      subtitle = str_wrap('Dots show true values, and lines involve interpolation for visualization purposes', width=70),
      caption = str_wrap(paste0('The ', unique(df_tall$age), ' yr-old ', 
                                unique(df_tall$sex), 
                                ' subject had a BMI of ', unique(round(df_tall$bmi, 1)), 
                                ' with waist and hip circumfurances (cm) of ', 
                                unique(round(df_tall$waist_circumf_cm, 1)), ' and ', 
                                unique(round(df_tall$hip_circumf_cm, 1)), ', respectively.') , width = 57)
    ) +
    
    ggplot2::theme_bw(base_size = 16) +
    ggplot2::theme(legend.position = 'none') +
    ggplot2::theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )+
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = 16),  # Increase subtitle size
      plot.caption = ggplot2::element_text(size = 16, hjust=0),  # Increase caption size
      axis.title = ggplot2::element_text(size = 14),  # Increase axis label size
      axis.text = ggplot2::element_text(size = 14) # Increase axis text size
    )

}



