{
    printf ("set pm3d map\n"); 
    # printf("set palette defined (-1 \"blue\", 0 \"white\", 1 \"red\")\n");
    printf("set cbrange [0:12000]\n");
    for(i=0.0;i<=20.0;i+=0.01){ printf("splot \"%010.2lf[orb].surfacedensity\" u 1:2:3 \n",i)}
}
