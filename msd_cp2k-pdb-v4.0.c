/*BY FILIPE MATUSALEM, SEPT 2020     fmatusa@unicamp.br */
/*Program to compute mean square displacement from CP2K PDB trajectory file*/
/*Compilation: g++ -o rmsd_cp2k-pdb.x rmsd_cp2k-pdb.c*/
//APRIL 2021  include coefficient diffusion calculation
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int main(int argc, char *argv[])
{
FILE *pdb,*meansquare,*meansquare2;
float a,b,c,timestep;
int i,j,k,l,m,natoms,nspecies,nsteps,ntype[10],range;
char str1[150],ch,species[10][10],lixo[150];


 if( argc < 3 ){
printf("\n\n");
printf("EEEEEEE   RRRRRRRR   RRRRRRRR   OOOOOOO   RRRRRRRR\n");
printf("EE        RR    RR   RR    RR   OO   OO   RR    RR\n");
printf("EE        RR    RR   RR    RR   OO   OO   RR    RR\n");
printf("EEEE      RRRRRRRR   RRRRRRRR   OO   OO   RRRRRRRR\n");
printf("EE        RRRR       RRRR       OO   OO   RRRR\n");
printf("EE        RR  RR     RR  RR     OO   OO   RR  RR\n");
printf("EEEEEEE   RR    RR   RR    RR   OOOOOOO   RR    RR\n\n");

printf("Enter the name of the cp2k trajectory pdb file and MD time step in fs as arguments\n\n");

 exit(0);}


printf("-----------------------------------------------------------------------------------------------\n\n");
printf("a third argument can be used to define the interval for multiple start points.\n");


strcpy(str1,argv[1]);

timestep=atof(argv[2]); 



pdb = fopen(str1,"r"); /* Arquivo ASCII, para leitura */
if(!pdb)
{
printf( "Error opening argument 1 file\n");
printf("Enter the name of the cp2k trajectory pdb file \n");

exit(0);
}

meansquare = fopen("meansquare.dat","w"); /* Arquivo ASCII, para escrita */
if(!meansquare)
{
printf( "Error creating meansquare.dat file\n");
exit(0);
}


do
fscanf(pdb,"%s",str1);                                      /*posiciona o  após a palavra CRYST1e*/
while(strcmp(str1,"CRYST1")!=0);

do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');

natoms=0;
do
{fscanf(pdb,"%s",str1);  if(strcmp(str1,"ATOM")==0)natoms++;  }                                   
while(strcmp(str1,"END")!=0);

printf("No atoms %d\n",natoms);
rewind(pdb);


do
fscanf(pdb,"%s",str1);                                      /*posiciona o  após a palavra CRYST1e*/
while(strcmp(str1,"CRYST1")!=0);

do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');

fscanf(pdb,"%s",str1);
fscanf(pdb,"%s",str1);
fscanf(pdb,"%s",species[0]);

do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');

j=k=1;
for(i=0;i<natoms-1;i++){
fscanf(pdb,"%s",str1);
fscanf(pdb,"%s",str1);
fscanf(pdb,"%s",species[j]);
k++;
if(strcmp(species[j-1],species[j])!=0){ntype[j-1]=k-1;k=1;j++;}

do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');
}
ntype[j-1]=k;
nspecies=j;

printf("No. Species %d \n\n",nspecies);
printf("Specie  Number\n");
for(i=0;i<nspecies;i++){
printf("   %s      %d \n",species[i],ntype[i]);
}
rewind(pdb);

nsteps=0;
while (fscanf(pdb,"%s",str1) != EOF){            /*conta steps*/
if(strcmp(str1,"CRYST1")==0)nsteps++;                      
}

printf("No. steps = %d\n",nsteps);
rewind(pdb);


float pos0[natoms][3],msd[natoms],msd_per_type[nspecies],msd_per_type_all[nspecies][nsteps];


do
fscanf(pdb,"%s",str1);                                      /*posiciona o  após a palavra CRYST1*/
while(strcmp(str1,"CRYST1")!=0);


do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');

for(i=0;i<natoms;i++){

fscanf(pdb,"%s",lixo);
fscanf(pdb,"%s",lixo);
fscanf(pdb,"%s",lixo);
fscanf(pdb,"%f",&pos0[i][0]);
fscanf(pdb,"%f",&pos0[i][1]);
fscanf(pdb,"%f",&pos0[i][2]);

do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');
}




fprintf(meansquare,"#  time (fs), msd (Ang^2) and diffusion D (cm^2 s^-1)   \n");
fprintf(meansquare,"#  time    ");
for(i=0;i<nspecies;i++)fprintf(meansquare,"    msd_%s       D_%s   ",species[i],species[i]);
fprintf(meansquare,"\n");


for(k=0;k<nsteps-1;k++){

for(i=0;i<natoms;i++)msd[i]=0;
for(i=0;i<nspecies;i++)msd_per_type[i]=0;


do
fscanf(pdb,"%s",str1);                                      /*posiciona o  após a palavra CRYST1*/
while(strcmp(str1,"CRYST1")!=0);

do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');

for(i=0;i<natoms;i++){

fscanf(pdb,"%s",lixo);
fscanf(pdb,"%s",lixo);
fscanf(pdb,"%s",lixo);
fscanf(pdb,"%f",&a); 
fscanf(pdb,"%f",&b); 
fscanf(pdb,"%f",&c); 

do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');

msd[i]=pow(sqrt(pow(pos0[i][0]-a,2)+pow(pos0[i][1]-b,2)+pow(pos0[i][2]-c,2)),2);
}

fprintf(meansquare,"  %f    ",(k+1)*timestep);

//average for each specie
i=0;
for(j=0;j<nspecies;j++){for(l=0;l<ntype[j];l++)msd_per_type[j]=msd_per_type[j]+msd[l+i]; 
i=i+l;

a=msd_per_type[j]/ntype[j]; msd_per_type_all[j][k]=a;
b=a*pow(10,-16);     //from Ang^2 to cm^2
b=b/(6*(k+1)*timestep*pow(10,-15));    //Diffusion
b=log10(b);   //Log_10 (Diffusion)


fprintf(meansquare,"%f  %f   ",msd_per_type_all[j][k],b);
}

fprintf(meansquare,"\n");

}








//---------------------------if the initial point varies --------
if(argc > 3){


range = atoi(argv[3]); printf("\nInterval for multiple start points: %d\n",range);

int nintervals;

nintervals = nsteps/range;
printf("\nNo. intervals: %d\n",nintervals);

float msd_per_type_all3d[nspecies][nsteps][nintervals], sum[nspecies][nsteps];

for(j=0;j<nspecies;j++)for(k=0;k<nsteps-1;k++)for(m=0;m<nintervals;m++)msd_per_type_all3d[j][k][m]=0;
for(j=0;j<nspecies;j++)for(k=0;k<nsteps-1;k++)sum[nspecies][nsteps]=0;

for(j=0;j<nspecies;j++)for(k=0;k<nsteps-1;k++)msd_per_type_all3d[j][k][0]=msd_per_type_all[j][k];

//loop para mudar o ponto de referencia inicial

for(m=1;m<nintervals;m++){

rewind(pdb);

for(i=0;i<m*range+1;i++){
do
fscanf(pdb,"%s",str1);                                      /*posiciona o  após a palavra CRYST1*/
while(strcmp(str1,"CRYST1")!=0);
}
//printf("Points excluded in interval %d: %d\n",m,i-1);



do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');

for(i=0;i<natoms;i++){

fscanf(pdb,"%s",lixo);
fscanf(pdb,"%s",lixo);
fscanf(pdb,"%s",lixo);
fscanf(pdb,"%f",&pos0[i][0]);
fscanf(pdb,"%f",&pos0[i][1]);
fscanf(pdb,"%f",&pos0[i][2]);

do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');
}





printf("Points included in interval %d: %d\n",m,nsteps-m*range-1);

for(k=0;k<nsteps-m*range-1;k++){

for(i=0;i<natoms;i++)msd[i]=0;
for(i=0;i<nspecies;i++)msd_per_type[i]=0;


do
fscanf(pdb,"%s",str1);                                      /*posiciona o  após a palavra CRYST1*/
while(strcmp(str1,"CRYST1")!=0);

do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');

for(i=0;i<natoms;i++){

fscanf(pdb,"%s",lixo);
fscanf(pdb,"%s",lixo);
fscanf(pdb,"%s",lixo);
fscanf(pdb,"%f",&a); 
fscanf(pdb,"%f",&b); 
fscanf(pdb,"%f",&c); 

do
ch = getc(pdb);              /*chega ao fim da linha*/
while(ch!='\n');

msd[i]=pow(sqrt(pow(pos0[i][0]-a,2)+pow(pos0[i][1]-b,2)+pow(pos0[i][2]-c,2)),2);
}



//average for each specie
i=0;
for(j=0;j<nspecies;j++){for(l=0;l<ntype[j];l++)msd_per_type[j]=msd_per_type[j]+msd[l+i]; 
i=i+l;

msd_per_type_all3d[j][k][m]=msd_per_type[j]/ntype[j]; 

}}
}

printf("\nMSD written to files:\n");

for(j=0;j<nspecies;j++)for(k=0;k<nsteps-1;k++)for(m=0;m<nintervals;m++)sum[j][k]=sum[j][k]+msd_per_type_all3d[j][k][m];


for(i=0;i<nspecies;i++){
sprintf(str1, "meansquare");
strcat(str1,"_");
strcat(str1,species[i]);
strcat(str1,".dat");
printf("%s\n",str1);

meansquare2 = fopen(str1,"w"); /* Arquivo ASCII, para escrita */
if(!meansquare2)
{
printf( "Error creating meansquare2.dat file\n");
exit(0);
}



fprintf(meansquare2,"#  time       ");
for(m=0;m<nintervals;m++)fprintf(meansquare2,"int_%d     ",m);
fprintf(meansquare2,"  average\n");





//for(k=0;k<nsteps-1;k++){fprintf(meansquare2,"  %f    ",(k+1)*timestep); for(m=0;m<nintervals;m++){fprintf(meansquare2,"%f  ",msd_per_type_all3d[i][k][m]);}fprintf(meansquare2,"\n");}

for(k=0;k<nsteps-1;k++){fprintf(meansquare2,"  %6.2f    ",(k+1)*timestep); for(m=0;m<nintervals;m++){fprintf(meansquare2,"%f  ",msd_per_type_all3d[i][k][m]);} l=(nsteps-2-k)/range+1;  fprintf(meansquare2,"  %f",sum[i][k]/l); fprintf(meansquare2,"\n");}
fclose(meansquare2);
}
}





fclose(pdb);
fclose(meansquare);

}
