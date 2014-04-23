#include <stdio.h>

int main(int argc, char *argv[])
{
  FILE *fp;
  fp = fopen("out.gp", "w");

  for(double s=0;s<200.0;s+=0.01){
    fprintf(fp, "set term post eps enhan color;set out \"%010.2lf.eps\"\n",s);
    fprintf(fp, "set pm3d map;splot \"%010.2lf[orb].surfacedensity\" u 1:2:3\n",s);
  }

  fclose(fp);
  return 0;
}
