/******************************************************************************/
/* Probability Weighted Iterative Generalized Least Squares                   */
/* for two level random coefficient models                                    */
/* necessary to declare:
/* string varlist :  covariates */
/* dep: dependent variable*/
/*  string varlist1: random effects */
/* cluster_var: cluster variable*/
/*  _wj_rep: cluster level weights*/
/* _wi_j : individual level weight*/
/******************************************************************************/

version 9

mata:

void pwigls_2l_adcv( string varlist,  dep, string varlist1, cluster_var,  _wj_rep,  _wi_j )
{
start= st_global("c(current_time)")
today= st_global("c(current_date)")
x = st_data(., tokens(varlist))
y = st_data(.,tokens(dep))
z = st_data(.,tokens(varlist1))
cluster_var = st_data(.,tokens(cluster_var))
cluster= uniqrows(cluster_var) 
wj_rep = st_data(.,tokens(_wj_rep))
wi_j = st_data(.,tokens(_wi_j))

nvar = cols(x)
ncluster = rows(cluster)
nsubjc_t=rows(y) 
q = cols(z)
s = ((q*(q+1))/2)+1

h_matrix= ((I(s)[vec(makesymmetric(invvech(1::s-1))), ]))'

name1 = tokens(varlist)
name3 = ("sigma2_u0","sigma_u10","sigma2_u1" ,"sigma2_e") 

lambidaj = J(nsubjc_t,1,0) 
info_j=panelsetup(cluster_var,1) 
cluster_wgt=J(ncluster,1,0)



/*-------------- Calculating the Scaled Weights --------------*/
/*-------------- Scaling Method 2               --------------*/

for (i=1; i<=rows(info_j) ; i++){
 nsubjc=panelsubmatrix(cluster_var, i, info_j)
 kwi_j=panelsubmatrix(wi_j, i, info_j)
 np = rows(nsubjc)
 k =mean(kwi_j)
 parte = J(np,1,k)

 a=info_j[i,1]
 b=info_j[i,2]
 lambidaj[a..b,]=parte 

 wj_rep_j=panelsubmatrix(wj_rep, i, info_j)

//cluster_wgt[i,]=mean(wj_rep_j)
cluster_wgt[i,]=uniqrows(wj_rep_j)

}
 
wi_j_star = wi_j :/ lambidaj //when colon is used it is element by element

wj_star = cluster_wgt / mean(cluster_wgt)
wjesc_r = wj_rep / mean(cluster_wgt)

/*-------------------------- Calculating Beta_zero ---------------------*/
mat_t1j = J(ncluster,nvar^2,0) // it is necessary to declare a matrix 
mat_t3j = J(ncluster,nvar,0)

for (i=1; i<=rows(info_j) ; i++){
 nsubjc=panelsubmatrix(cluster_var, i, info_j)
 yj=panelsubmatrix(y, i, info_j)
 xj=panelsubmatrix(x, i, info_j)
 np = rows(nsubjc)
 
 wi_j_starj=panelsubmatrix(wi_j_star, i, info_j)
 diag = diag(wi_j_starj)
 wj = wj_star[i,]

// each cluster is one line of the matrices 
 mat_t1j[i,] = (vec( ( xj' * diag * xj ) * wj ))'
 mat_t3j[i,] = (vec( ( xj' * diag * yj ) * wj))'
}

// Calculates the partial sums
somat1 = rowshape(colsum(mat_t1j),nvar)
somat3 = colsum(mat_t3j)

beta0 = cholinv(somat1) * (somat3')


/*-------Calculating - v(u0j) e v(eij)--------------- */
vec_t6 = J(ncluster,1,0)
vec_aux = J(ncluster,1,0)
teta0 =  (vech(diag(0.5):*I(q))\0)

for (i=1; i<=rows(info_j) ; i++){
 nsubjc=panelsubmatrix(cluster_var, i, info_j)
 yj=panelsubmatrix(y, i, info_j)
 xj=panelsubmatrix(x, i, info_j)
 np = rows(nsubjc)

 wi_j_starj=panelsubmatrix(wi_j_star, i, info_j)
 wj = wj_star[i,]

//raw residuals
 resid = yj - xj * beta0

//cluster 
 uj0 = (wi_j_starj' * resid ) / colsum(wi_j_starj)

 vec_t6[i,] = (wi_j_starj' * ( (resid :- uj0) :^ 2 ) ) * wj
 vec_aux[i,]= wj * (colsum(wi_j_starj)- 1 )
}
 teta0[s,]=colsum(vec_t6)/colsum(vec_aux)



/*--------------------------------------------------------------------------
                      PWIGLS - ITERATIVE MODULE
--------------------------------------------------------------------------*/

matp = J(ncluster,nvar^2,0)
matq = J(ncluster,nvar,0)

beta_ant = beta0
beta = beta_ant:*2 //It is a trick

teta_ant=teta0
teta=teta_ant :* 2 //it is a trick 


itera = 1
 

// It will do up to 200 iterations where precision is 0.000001

while (itera<=200 & (any(abs((teta-teta_ant)):> 0.000001) | any(abs((beta-beta_ant)):> 0.000001) ))
{



 /*--------------looping for beta---------------------------*/
for (i=1; i<=rows(info_j) ; i++){
 nsubjc=panelsubmatrix(cluster_var, i, info_j)
 yj=panelsubmatrix(y, i, info_j)
 xj=panelsubmatrix(x, i, info_j)
 zj=panelsubmatrix(z, i, info_j)
 
 np = rows(nsubjc)
 
 v=panelsubmatrix(wi_j_star, i, info_j)
 diag = diag(v)

 wj = wj_star[i,]

 t1j = ( xj' * diag * xj )
 t2j = ( xj' * diag * zj )
 t3j = ( xj' * diag * yj )
 t4j = ( zj' * diag * yj )
 t5j = ( zj' * diag * zj )

 if (itera ~= 1) {
 aj= cholinv( t5j + teta[s,]:*(cholinv(invvech(teta[1..(s-1),]))))
 }
 else {
 aj= cholinv( t5j + teta0[s,]:*(cholinv(invvech(teta0[1..(s-1),]))))
 }

 matp[i,] = (vec ( wj :* (t1j - t2j * aj * t2j') ) )'
 matq[i,] = (vec ( wj :* (t3j - t2j * aj * t4j ) ) )'
}

/*----------- beta-------------- */
s_matp = rowshape(colsum(matp),nvar)
s_matq = colsum(matq)

if (itera ~= 1){
beta_ant = beta
 }

beta = cholinv(s_matp) * (s_matq')



/*--------- looping for teta=inv(R)x S ------------------*/
 q = cols(z)
 s = ((q*(q+1))/2)+1

R = J(ncluster,s*s,0)
S = J(ncluster,s,0)
 

for (i=1; i<=rows(info_j) ; i++){
 nsubjc=panelsubmatrix(cluster_var, i, info_j)
 yj=panelsubmatrix(y, i, info_j)
 xj=panelsubmatrix(x, i, info_j)
 zj=panelsubmatrix(z, i, info_j)

 np = rows(nsubjc)

 v=panelsubmatrix(wi_j_star, i, info_j)
 diag = diag(v)
 
 wi_j_starj=panelsubmatrix(wi_j_star, i, info_j)
 wj = wj_star[i,]

 t1j = ( xj' * diag * xj )
 t2j = ( xj' * diag * zj )
 t3j = ( xj' * diag * yj )
 t4j = ( zj' * diag * yj )
 t5j = ( zj' * diag * zj )
 eij = yj - xj * beta

Rklj = J(s,s,0)
Skj = J(s,1,0)
B= J(s,q*q,0)
C = J(s,q*q,0)
//H = J(s,q*q,0)


H = h_matrix


delta = J(1,s,0)
delta[.,s]=1
 

 if (itera ~= 1) 
{
aj= cholinv( t5j + teta[s,]:*(cholinv(invvech(teta[1..(s-1),]))))

for (k=1; k<=s ; k++){

B[k,]=(vec(teta[s,]*aj*cholinv(invvech(teta[1..(s-1),]))* rowshape(H[k,],q)- delta[,k]*aj ))'

C[k,]=(vec(-delta[,k]*aj + rowshape(B[k,],q)' - rowshape(B[k,],q)'* t5j * aj ))'

for (l=1; l<=s ; l++){

Rklj[k,l]= wj*(delta[.,k]*delta[.,l]*colsum(v) + delta[.,l]*trace(t5j*rowshape(C[k,],q)')+ delta[.,k]*trace(t5j*rowshape(H[l,],q))+trace(t5j*rowshape(C[k,],q)'*t5j*rowshape(H[l,],q)))

} //end of l

Skj[k,]= wj*(delta[.,k]*trace(eij'*diag*eij) +trace(eij'*diag*zj * rowshape(C[k,],q)' * zj'*diag*eij))
 } //end of k
 } //end of if
 else 
{
 aj= cholinv( t5j + teta0[s,]:*(cholinv(invvech(teta0[1..(s-1),]))))

for (k=1; k<=s ; k++){
B[k,]=(vec( teta0[s,]*aj*cholinv(invvech(teta0[1..(s-1),]))* rowshape(H[k,],q)- delta[,k]*aj ))'
C[k,]=(vec(-delta[,k]*aj + rowshape(B[k,],q)' - rowshape(B[k,],q)'*t5j*aj))'
for (l=1; l<=s ; l++){

Rklj[k,l]= wj*(delta[.,k]*delta[.,l]*colsum(v) + delta[.,l]*trace(t5j*rowshape(C[k,],q)')+ delta[.,k]*trace(t5j*rowshape(H[l,],q))+trace(t5j*rowshape(C[k,],q)'*t5j*rowshape(H[l,],q)))

} //end of l
Skj[k,]= wj*(delta[.,k]*trace(eij'*diag*eij) +trace(eij'*diag*zj * rowshape(C[k,],q)' * zj'*diag*eij))
} //end of k
} //end of else

R[i,]=(vec(Rklj))'
S[i,]=(vec(Skj))'
} //end of for

 matr=colsum(R)
 mats=colsum(S)

 r_mat = rowshape(matr,s)
 s_mat = rowshape(mats,s)

  
/*-----teta = u0j , eij----*/

if (itera ~= 1) {
 teta_ant = teta
		}

teta = cholinv(r_mat) * s_mat

 

/*------End of iterative process----------*/
itera = itera + 1
}

/* Number of Iterations*/
n_it = itera - 1
 
  

/*-------------------------------------------------------------------------
 Residuals
-------------------------------------------------------------------------*/
//level 2
u = J(ncluster,q, 0)  
var_u = J(ncluster,q*q,0)  
dp_u = J(ncluster,q,0)  



yhat = J(nsubjc_t,1, 0)

//level1
v = J(nsubjc_t,1, 0)

var_v = J(nsubjc_t,1, 0)


for (i=1; i<=rows(info_j) ; i++){

 nsubjc=panelsubmatrix(cluster_var, i, info_j)
 yj=panelsubmatrix(y, i, info_j)
 xj=panelsubmatrix(x, i, info_j)
 zj=panelsubmatrix(z, i, info_j)

 np = rows(nsubjc)

 ej = yj - xj * beta

 a=info_j[i,1]
 b=info_j[i,2]
 
 Sigmau=(invvech(teta[1..(s-1),]))

 Rhj = Sigmau* zj'

 Vj = zj * Sigmau * zj' + teta[s,] :* I(np)

 aux =  Rhj * cholinv(Vj)
 

 u[i,] = (aux * ej)'
 
//In goldstein we multiply this by diag if this is different than the gllamm one this is the reason why

 var_u[i,] = (vec(Sigmau - aux * Rhj'))'
 dp_u[i,] = (vec(sqrt(diagonal(rowshape(var_u[i,],q) ))))'

// it is also different than goldstein 

 yhat1 = xj * beta + zj*u[i,]'
 yhat[a..b,] = yhat1

 vj = ej - zj*u[i,]'
 v[a..b,] = vj

 aux = teta[s,] :* ( 1 - (1/np) )

 var = J(np,1,aux)
 var_v[a..b,] = var 
}

u_pad = u :/  dp_u

v_pad = v :/ sqrt(var_v)
 



/*-------------------------------------------------------------------------
 Variances of Beta and Teta
-------------------------------------------------------------------------*/
mat_c = J(ncluster,nvar^2,0)
mat_d = J(ncluster,s*s,0)

for (i=1; i<=rows(info_j) ; i++){
 nsubjc=panelsubmatrix(cluster_var, i, info_j)
 yj=panelsubmatrix(y, i, info_j)
 xj=panelsubmatrix(x, i, info_j)
 zj=panelsubmatrix(z, i, info_j)

 np = rows(nsubjc)
 wj = wj_star[i,]
 v = panelsubmatrix(wi_j_star, i, info_j)
 diag = diag(v)

/*----------Beta----------*/

 ej = yj - xj * beta

 t2j = ( xj' * diag * zj )
 t4j = ( zj' * diag * yj )
 t5j = ( zj' * diag * zj )
 t7j = ( xj' * diag * ej )
 t8j = ( ej' * diag * zj )
 
 aj= cholinv( t5j + teta[s,]*(cholinv(invvech(teta[1..(s-1),]))))
 
// cj = t7j - ( aj * t2j * t8j)
 cj = t7j - ( t2j * aj *  t8j')

 mat_c[i,] = (vec( (wj:^2) * (cj * cj') ))' 
 
/*----------teta----------*/

 Rklj=rowshape(R[i,],s)
 Skj=rowshape(S[i,],s)

 Dkj =  (Skj-Rklj*teta)*(Skj-Rklj*teta)'
 mat_d[i,]=(vec(Dkj))'
}
 

/*----------Variances-----------------*/
var_beta = cholinv(s_matp)*((ncluster /( ncluster-1))*rowshape(colsum(mat_c),nvar))*cholinv(s_matp)
dp_beta = sqrt(diagonal(var_beta))
//dp_beta

var_teta = cholinv(r_mat)*ncluster/(ncluster-1)* rowshape(colsum(mat_d),s)* cholinv(r_mat)
dp_teta = sqrt(diagonal(var_teta))

finish= st_global("c(current_time)")
today1= st_global("c(current_date)")

z_star = beta:/dp_beta
z_star2= teta:/dp_teta

z_star_l= beta-abs(invnormal(0.025)):*dp_beta
z_star_u= beta+abs(invnormal(0.025)):*dp_beta
pz_star= 2:*(1:-normal(abs(z_star)))

z_star2_l= teta-abs(invnormal(0.05)):*dp_teta
z_star2_u= teta+abs(invnormal(0.05)):*dp_teta
pz_star2= 2:*(1:-normal(abs(z_star2)))

st_matrix("beta",beta)
st_matrix("var_beta",var_beta)
 
//result1 = (beta , dp_beta , z_star)
//result2 = (teta , dp_teta , z_star2)
//name2
//name1' 
//result1
//name3'
//result2

//sqrt(diagonal(cholinv(1/teta[s,]:*s_matp)))
//sqrt(diagonal(2*cholinv(1/(teta[s,]^2):*r_mat)))


printf("\n")
printf("{hline 87}\n")
printf("            Probability Weighted Iterative Generalized Least Squares            \n")
printf("{hline 87}\n")
printf("General Information\n")
printf("\n")
printf("Response Variable = {txt}%19s \n", dep)
printf("Weight at Level 2 = {txt}%19s \n", _wj_rep)
printf("Weight at Level 1 = {txt}%19s \n", _wi_j)
printf("\n")
printf("Start running on %s at  %s\n",  today , start)
printf("Number of Iterations    = %3.0f\n", n_it)
printf("Number of Level 1 units = %3.0f\n", nsubjc_t)
printf("Number of Level 2 units = %3.0f\n", ncluster)
printf("\n")
printf("{hline 87}\n")
printf("       Fixed Effects{c |}   Coef.    Std.Err.    z     P>|z|   [95%sConf.Interval]  Init.Val.\n","%")
printf("{hline 20}{c +}{hline 66}\n")
for (mi=1; mi<= nvar; mi++) {
printf("{txt}%19s {c |} {res}%8.0g  %8.0g  %6.2f  %6.3f  %8.0g  %8.0g  %8.0g\n",name1[.,mi]', beta[mi,.], dp_beta[mi,.],z_star[mi,.], pz_star[mi,.] , z_star_l[mi,.],z_star_u[mi,.],beta0[mi,.] )
}
printf("{hline 87}\n")
printf("\n")

printf("{hline 87}\n")
printf(" Variance Components{c |}   Coef.    Std.Err.    z     P>|z|   [95%sConf.Interval]  Init.Val.\n","%")
printf("{hline 20}{c +}{hline 66}\n")
for (mi=1; mi<= s; mi++) {
printf("{txt}%19s {c |} {res}%8.0g  %8.0g  %6.2f  %6.3f  %8.0g  %8.0g  %8.0g\n",name3[.,mi]', teta[mi,.], dp_teta[mi,.],z_star2[mi,.], pz_star2[mi,.] , z_star2_l[mi,.],z_star2_u[mi,.],teta0[mi,.] )

}
printf("{hline 87}\n")
printf("Note: Robust Standard Errors \n")
printf("      Finish running on %s at  %s\n",  today1 , finish)
printf("{hline 87}\n")

}

//if (direxists("c:\ado")==0) _mkdir("c:\ado")
//if (direxists("c:\ado\PERSONAL")==0) _mkdir("c:\ado\PERSONAL")
bb2=0		/* dummy command to end the if structure */
mata mosave  pwigls_2l_adcv(), replace

end
