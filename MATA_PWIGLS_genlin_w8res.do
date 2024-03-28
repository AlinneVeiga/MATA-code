/********************************************************************************/
/* Probability Weighted Iterative Generalized Least Squares                   	*/
/* for multivariate multilevel models - for longitudinal data                 	*/
/* necessary to declare:							*/
/* string varlist :  covariates 						*/
/* dep: dependent variable							*/
/* string varlist1: random effects 						*/
/* string varlist2: time variables 						*/
/* cluster_var: cluster variable						*/
/* _wj_rep: cluster level weights						*/
/* _wi_j : individual level weight						*/
/******************************************************************************/
version 9

mata:

void pwigls_genlin_adcv( string varlist,  dep, string varlist2, string varlist1, cluster_var,  _wj_rep,  _wi_j )

{
start= st_global("c(current_time)")
today= st_global("c(current_date)")
x = st_data(., tokens(varlist))
y = st_data(.,tokens(dep))
z = st_data(.,tokens(varlist1))
time_var = st_data(.,tokens(varlist2))
cluster_var = st_data(.,tokens(cluster_var))
cluster= uniqrows(cluster_var) 
wj_rep = st_data(.,tokens(_wj_rep))
wi_j = st_data(.,tokens(_wi_j))


//the different numbers needed
nvar = cols(x)
ncluster = rows(cluster)
nsubjc_t=rows(y) // this will need change for when 3 levels or longitudinal

q = cols(z)

t = cols(time_var)

s = 12

h_matrix= ((I(s)[1, ]))'
h_matrix

delta_matrix=J(s,64,0)

delta_matrix[2,]=(1,0,0,0,0,0,0,0,
0,1,0,0,0,0,0,0,
0,0,1,0,0,0,0,0,
0,0,0,1,0,0,0,0,
0,0,0,0,1,0,0,0,
0,0,0,0,0,1,0,0,
0,0,0,0,0,0,1,0,
0,0,0,0,0,0,0,1)


delta_matrix[3,]=(1,1,0,0,0,0,0,0,
1,1,1,0,0,0,0,0,
0,1,1,1,0,0,0,0,
0,0,1,1,0,0,0,0,
0,0,0,0,1,1,0,0,
0,0,0,0,1,1,1,0,
0,0,0,0,0,1,1,1,
0,0,0,0,0,0,1,1)

delta_matrix[4,]=(1,0,1,0,0,0,0,0,
0,1,0,1,0,0,0,0,
1,0,1,0,0,0,0,0,
0,1,0,1,0,0,0,0,
0,0,0,0,1,0,1,0,
0,0,0,0,0,1,0,1,
0,0,0,0,1,0,1,0,
0,0,0,0,0,1,0,1)

delta_matrix[5,]=(1,0,0,1,0,0,0,0,
0,1,0,0,0,0,0,0,
0,0,1,0,0,0,0,0,
1,0,0,1,0,0,0,0,
0,0,0,0,1,0,0,1,
0,0,0,0,0,1,0,0,
0,0,0,0,0,0,1,0,
0,0,0,0,1,0,0,1)

delta_matrix[6,]=(1,0,0,0,0,0,0,0,
0,1,0,0,0,0,0,0,
0,0,1,0,0,0,0,0,
0,0,0,1,1,0,0,0,
0,0,0,1,1,0,0,0,
0,0,0,0,0,1,0,0,
0,0,0,0,0,0,1,0,
0,0,0,0,0,0,0,1)

delta_matrix[7,]=(1,0,0,0,0,0,0,0,
0,1,0,0,0,0,0,0,
0,0,1,0,1,0,0,0,
0,0,0,1,0,1,0,0,
0,0,1,0,1,0,0,0,
0,0,0,1,0,1,0,0,
0,0,0,0,0,0,1,0,
0,0,0,0,0,0,0,1)

delta_matrix[8,]=(1,0,0,0,0,0,0,0,
0,1,0,0,1,0,0,0,
0,0,1,0,0,1,0,0,
0,0,0,1,0,0,1,0,
0,1,0,0,1,0,0,0,
0,0,1,0,0,1,0,0,
0,0,0,1,0,0,1,0,
0,0,0,0,0,0,0,1)
 
delta_matrix[9,]=(1,0,0,0,1,0,0,0,
0,1,0,0,0,1,0,0,
0,0,1,0,0,0,1,0,
0,0,0,1,0,0,0,1,
1,0,0,0,1,0,0,0,
0,1,0,0,0,1,0,0,
0,0,1,0,0,0,1,0,
0,0,0,1,0,0,0,1)

delta_matrix[10,]=(1,0,0,0,0,1,0,0,
0,1,0,0,0,0,1,0,
0,0,1,0,0,0,0,1,
0,0,0,1,0,0,0,0,
0,0,0,0,1,0,0,0,
1,0,0,0,0,1,0,0,
0,1,0,0,0,0,1,0,
0,0,1,0,0,0,0,1)

delta_matrix[11,]=(1,0,0,0,0,0,1,0,
0,1,0,0,0,0,0,1,
0,0,1,0,0,0,0,0,
0,0,0,1,0,0,0,0,
0,0,0,0,1,0,0,0,
0,0,0,0,0,1,0,0,
1,0,0,0,0,0,1,0,
0,1,0,0,0,0,0,1)

delta_matrix[12,]=(1,0,0,0,0,0,0,1,
0,1,0,0,0,0,0,0,
0,0,1,0,0,0,0,0,
0,0,0,1,0,0,0,0,
0,0,0,0,1,0,0,0,
0,0,0,0,0,1,0,0,
0,0,0,0,0,0,1,0,
1,0,0,0,0,0,0,1)

rowshape(delta_matrix[1,],8)
rowshape(delta_matrix[2,],8)
rowshape(delta_matrix[3,],8)
rowshape(delta_matrix[4,],8)
rowshape(delta_matrix[5,],8)
rowshape(delta_matrix[6,],8)
rowshape(delta_matrix[7,],8)
rowshape(delta_matrix[8,],8)
rowshape(delta_matrix[9,],8)
rowshape(delta_matrix[10,],8)
rowshape(delta_matrix[11,],8)
rowshape(delta_matrix[12,],8)


name1 = tokens(varlist)
name3 = ("Sigma_u_2","Genlin(1)","Genlin(2)","Genlin(3)","Genlin(4)","Genlin(5)","Genlin(6)","Genlin(7)","Genlin(8)","Genlin(9)","Genlin(10)","Genlin(11)","Genlin(12)") 
 

lambidaj = J(nsubjc_t,1,0) // Lambda_j to be used in the scaling method

info_j=panelsetup(cluster_var,1) //This is to be used when doing procedures for each of the clusters

cluster_wgt=J(ncluster,1,0)

/*-------------- Calculating the Scaled Weights --------------*/
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

beta0 = invsym(somat1) * (somat3')

/*-------Calculating - v(u0j) e v(eij)--------------- */
vec_t6 = J(ncluster,1,0)
vec_aux = J(ncluster,1,0)
teta0 =  J(s,1,.5) 


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
 teta0[2,]=colsum(vec_t6)/colsum(vec_aux)


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

 if (itera ~= 1) {
 teta_genlin=teta[2,]*rowshape(delta_matrix[2,],8)+ teta[3,]*rowshape(delta_matrix[3,],8)+ teta[4,]*rowshape(delta_matrix[4,],8)+ teta[5,]*rowshape(delta_matrix[5,],8)+ teta[6,]*rowshape(delta_matrix[6,],8)+teta[7,]*rowshape(delta_matrix[7,],8)+teta[8,]*rowshape(delta_matrix[8,],8)+teta[9,]*rowshape(delta_matrix[9,],8)+teta[10,]*rowshape(delta_matrix[10,],8)+teta[11,]*rowshape(delta_matrix[11,],8)+teta[12,]*rowshape(delta_matrix[12,],8) 
 sigma=I(np/t)#teta_genlin
 aj= invsym(invsym(teta[1,])+zj'*(diag*invsym(sigma))*zj)
 invvj=diag*invsym(sigma)-diag*invsym(sigma)*zj*aj*zj'*diag*invsym(sigma)
 }
 else {
 teta0_genlin=teta0[2,]*rowshape(delta_matrix[2,],8)+ teta0[3,]*rowshape(delta_matrix[3,],8)+ teta0[4,]*rowshape(delta_matrix[4,],8)+ teta0[5,]*rowshape(delta_matrix[5,],8)+ teta0[6,]*rowshape(delta_matrix[6,],8)+teta0[7,]*rowshape(delta_matrix[7,],8)+teta0[8,]*rowshape(delta_matrix[8,],8)+teta0[9,]*rowshape(delta_matrix[9,],8)+teta0[10,]*rowshape(delta_matrix[10,],8)+teta0[11,]*rowshape(delta_matrix[11,],8)+teta0[12,]*rowshape(delta_matrix[12,],8)
 sigma0=I(np/t)#teta0_genlin
 aj=invsym(invsym(teta0[1,])+zj'*(diag*invsym(sigma0))*zj)
 invvj=diag*invsym(sigma0)-diag*invsym(sigma0)*zj*aj*zj'*diag*invsym(sigma0)
 }

 matp[i,] = (vec ( wj :* xj'*invvj*xj ))'
 matq[i,] = (vec ( wj :* xj'*invvj*yj ))'

}


/*----------- beta-------------- */
s_matp = rowshape(colsum(matp),nvar)
s_matq = colsum(matq)

if (itera ~= 1){
beta_ant = beta
 }

beta = invsym(s_matp) * (s_matq')



/*--------- looping for teta=inv(R)x S ------------------*/
 q = cols(z)
 s = 12
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
 wj = wj_star[i,]
 ej = yj - xj * beta

Rklj = J(s,s,0)
Skj = J(s,1,0)

 H = h_matrix
 delta=delta_matrix


 if (itera ~= 1) {
 teta_genlin=teta[2,]*rowshape(delta_matrix[2,],8)+ teta[3,]*rowshape(delta_matrix[3,],8)+ teta[4,]*rowshape(delta_matrix[4,],8)+ teta[5,]*rowshape(delta_matrix[5,],8)+ teta[6,]*rowshape(delta_matrix[6,],8)+teta[7,]*rowshape(delta_matrix[7,],8)+teta[8,]*rowshape(delta_matrix[8,],8)+teta[9,]*rowshape(delta_matrix[9,],8)+teta[10,]*rowshape(delta_matrix[10,],8)+teta[11,]*rowshape(delta_matrix[11,],8)+teta[12,]*rowshape(delta_matrix[12,],8) 
 sigma=I(np/t)#teta_genlin
 aj= invsym(invsym(teta[1,])+zj'*(diag*invsym(sigma))*zj)
 invvj=diag*invsym(sigma)-diag*invsym(sigma)*zj*aj*zj'*diag*invsym(sigma)
 }
 else {
 teta0_genlin=teta0[2,]*rowshape(delta_matrix[2,],8)+ teta0[3,]*rowshape(delta_matrix[3,],8)+ teta0[4,]*rowshape(delta_matrix[4,],8)+ teta0[5,]*rowshape(delta_matrix[5,],8)+ teta0[6,]*rowshape(delta_matrix[6,],8)+teta0[7,]*rowshape(delta_matrix[7,],8)+teta0[8,]*rowshape(delta_matrix[8,],8)+teta0[9,]*rowshape(delta_matrix[9,],8)+teta0[10,]*rowshape(delta_matrix[10,],8)+teta0[11,]*rowshape(delta_matrix[11,],8)+teta0[12,]*rowshape(delta_matrix[12,],8) 
 sigma0=I(np/t)#teta0_genlin
 aj=invsym(invsym(teta0[1,])+zj'*(diag*invsym(sigma0))*zj)
 invvj=diag*invsym(sigma0)-diag*invsym(sigma0)*zj*aj*zj'*diag*invsym(sigma0)
 }
 

for (k=1; k<=s ; k++)
{
for (l=1; l<=s ; l++)
{
Rklj[k,l]= wj*(trace((invvj*(zj*h_matrix[k,]*zj' + I(np/t)#rowshape(delta_matrix[k,],t)' * invsym(diag))) * (invvj*(zj*h_matrix[l,]*zj' + I(np/t)#rowshape(delta_matrix[l,],t)'*invsym(diag)))))
} //end of l
Skj[k,]= wj*(trace(invvj*(zj*h_matrix[k,]*zj' + I(np/t)#rowshape(delta_matrix[k,],t)'*invsym(diag))*invvj*(ej*ej')))
 }   

R[i,]=(vec(Rklj))'
S[i,]=(vec(Skj))'
} 

 matr=colsum(R)
 mats=colsum(S)

 r_mat = rowshape(matr,s)
 s_mat = rowshape(mats,s)

  
/*-----teta = u0j , eij----*/

if (itera ~= 1) {
 teta_ant = teta
		}

teta = invsym(r_mat) * s_mat
 
teta
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
 
 Sigmau=(invvech(teta[1,]))

 Rhj = Sigmau* zj'

 teta_genlin=teta[2,]*rowshape(delta_matrix[2,],8)+ teta[3,]*rowshape(delta_matrix[3,],8)+ teta[4,]*rowshape(delta_matrix[4,],8)+ teta[5,]*rowshape(delta_matrix[5,],8)+ teta[6,]*rowshape(delta_matrix[6,],8)+teta[7,]*rowshape(delta_matrix[7,],8)+teta[8,]*rowshape(delta_matrix[8,],8)+teta[9,]*rowshape(delta_matrix[9,],8)+teta[10,]*rowshape(delta_matrix[10,],8)+teta[11,]*rowshape(delta_matrix[11,],8)+teta[12,]*rowshape(delta_matrix[12,],8) 
 sigma=I(np/t)#teta_genlin

 Vj = zj * Sigmau * zj' + sigma

 aux =  Rhj * invsym(Vj)
 
 u[i,] = (aux * ej)'
 
//In goldstein we multiply this by diag if this is different than the gllamm one this is the reason why

 var_u[i,] = (vec(Sigmau - aux * Rhj'))'
 dp_u[i,] = (vec(sqrt(diagonal(rowshape(var_u[i,],q) ))))'

// it is also different than goldstein 

 yhat1 = xj * beta + zj*u[i,]'
 yhat[a..b,] = yhat1

 vj = ej - zj*u[i,]'
 v[a..b,] = vj

// aux = sqrt(diagonal(I(np/t)#Toeplitz(teta[2..s,],teta[2..s,]'):* ( 1 - (1/(np/t)) )))
// var = J(np,1,aux)
//var_v[a..b,] = aux 
}

//u_pad = u :/  dp_u
//v_pad = v :/(var_v) //jah eh sqrt

//st_matrix("u",u)

resindex01 = st_addvar("float","u")
st_store((1,rows(u)),resindex01,u)


resindex0 = st_addvar("float","se")
st_store((1,rows(dp_u)),resindex0,dp_u)

resindex = st_addvar("float","resid")
st_store((1,rows(res)),resindex,res)

resindex1 = st_addvar("float","yhat1")
st_store((1,rows(yhat)),resindex1,yhat)

resindex2 = st_addvar("float","clusteru")
st_store((1,rows(cluster)),resindex2,cluster)

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

 teta_genlin=teta[2,]*rowshape(delta_matrix[2,],8)+ teta[3,]*rowshape(delta_matrix[3,],8)+ teta[4,]*rowshape(delta_matrix[4,],8)+ teta[5,]*rowshape(delta_matrix[5,],8)+ teta[6,]*rowshape(delta_matrix[6,],8)+teta[7,]*rowshape(delta_matrix[7,],8)+teta[8,]*rowshape(delta_matrix[8,],8)+teta[9,]*rowshape(delta_matrix[9,],8)+teta[10,]*rowshape(delta_matrix[10,],8)+teta[11,]*rowshape(delta_matrix[11,],8)+teta[12,]*rowshape(delta_matrix[12,],8) 
 sigma=I(np/t)#teta_genlin
 aj= invsym(invsym(teta[1,])+zj'*(diag*invsym(sigma))*zj)
 invvj=diag*invsym(sigma)-diag*invsym(sigma)*zj*aj*zj'*diag*invsym(sigma)
 
  
 cj = xj'*invvj*ej 

 mat_c[i,] = (vec( (wj:^2) * (cj * cj') ))' 

/*----------teta----------*/

 Rklj=rowshape(R[i,],s)
 Skj=rowshape(S[i,],s)

 Dkj =  (Skj-Rklj*teta)*(Skj-Rklj*teta)'
 mat_d[i,]=(vec(Dkj))'
}
 

/*----------Variances-----------------*/
var_beta = invsym(s_matp)*((ncluster /( ncluster-1))*rowshape(colsum(mat_c),nvar))*invsym(s_matp)
dp_beta = sqrt(diagonal(var_beta))


var_teta = invsym(r_mat)*ncluster/(ncluster-1)* rowshape(colsum(mat_d),s)* invsym(r_mat)
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

teta_genlin=teta[2,]*rowshape(delta_matrix[2,],8)+ teta[3,]*rowshape(delta_matrix[3,],8)+ teta[4,]*rowshape(delta_matrix[4,],8)+ teta[5,]*rowshape(delta_matrix[5,],8)+ teta[6,]*rowshape(delta_matrix[6,],8)+teta[7,]*rowshape(delta_matrix[7,],8)+teta[8,]*rowshape(delta_matrix[8,],8)+teta[9,]*rowshape(delta_matrix[9,],8)+teta[10,]*rowshape(delta_matrix[10,],8)+teta[11,]*rowshape(delta_matrix[11,],8)+teta[12,]*rowshape(delta_matrix[12,],8) 
TOEP=teta_genlin
Sigma_r=teta[1]*J(t,t,1)+ TOEP

st_matrix("beta",beta)
st_matrix("var_beta",var_beta)


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
printf("Number of Time points   = %3.0f\n", t)
printf("Number of Level 1 units = %3.0f\n", nsubjc_t/t)
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
printf("\n\n General Linear Matrix\n")
TOEP
printf("\n\n Total Variance\n")
Sigma_r
printf("{hline 87}\n")
printf("Note: Robust Standard Errors \n")
printf("      Finish running on %s at  %s\n",  today1 , finish)
printf("{hline 87}\n")

}

//if (direxists("c:\ado")==0) _mkdir("c:\ado")
//if (direxists("c:\ado\PERSONAL")==0) _mkdir("c:\ado\PERSONAL")
bb2=0		/* dummy command to end the if structure */
mata mosave  pwigls_genlin_adcv(),  replace

end
 
