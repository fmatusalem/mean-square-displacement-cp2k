/*BY FILIPE MATUSALEM, SEPT 2020     filipematus@gmail.com */
/*Program to compute root mean square displacement from CP2K PDB trajectory file*/
/*Compilation: g++ -o rmsd_cp2k-pdb.x rmsd_cp2k-pdb.c*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


int main(int argc, char *argv[])
{
FILE *pdb,*meansquare;
float a,b,c;
int i,j,k,l,natoms,nspecies,nsteps,ntype[10];
char str1[150],ch,species[10][10],lixo[150];


 if( argc < 2 ){
printf("\n\n");
printf("EEEEEEE   RRRRRRRR   RRRRRRRR   OOOOOOO   RRRRRRRR\n");
printf("EE        RR    RR   RR    RR   OO   OO   RR    RR\n");
printf("EE        RR    RR   RR    RR   OO   OO   RR    RR\n");
printf("EEEE      RRRRRRRR   RRRRRRRR   OO   OO   RRRRRRRR\n");
printf("EE        RRRR       RRRR       OO   OO   RRRR\n");
printf("EE        RR  RR     RR  RR     OO   OO   RR  RR\n");
printf("EEEEEEE   RR    RR   RR    RR   OOOOOOO   RR    RR\n\n");

printf("Enter the name of the cp2k trajectory pdb file \n\n");

 exit(0);}


printf("-----------------------------------------------------------------------------------------------\n\n");



strcpy(str1,argv[1]);

pdb = fopen(str1,"r"); /* Arquivo ASCII, para leitura */
if(!pdb)
{
printf( "Error opening argument 1 file\n");
printf("Enter the name of the cp2k trajectory pdb file (without .pdb)\n");

exit(0);
}

meansquare = fopen("rootmeansquare.dat","w"); /* Arquivo ASCII, para escrita */
if(!meansquare)
{
printf( "Error creating rootmeansquare.dat file\n");
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

float pos0[natoms][3],msd[natoms],msd_per_type[nspecies];

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





fprintf(meansquare,"#Step   ");for(i=0;i<nspecies;i++)fprintf(meansquare,"  %s       ",species[i]);fprintf(meansquare,"\n");

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

msd[i]=msd[i]+pow(sqrt(pow(pos0[i][0]-a,2)+pow(pos0[i][1]-b,2)+pow(pos0[i][2]-c,2)),2);
}

fprintf(meansquare,"  %d    ",k+1);

//average for each specie
i=0;
for(j=0;j<nspecies;j++){for(l=0;l<ntype[j];l++)msd_per_type[j]=msd_per_type[j]+msd[l+i];
i=i+l;

fprintf(meansquare,"%f  ",sqrt(msd_per_type[j]/ntype[j]));
}

fprintf(meansquare,"\n");

}



fclose(pdb);
fclose(meansquare);
}
