#include <stdio.h>
#include <stdlib.h>
#include <xdrfile/xdrfile_xtc.h>
#include <math.h>

int main ( int argc, char *argv[])
{
    if ( argc != 2 ) /* argc should be 2 for correct execution */
        {
		printf( "usage: %s filename\n", argv[0] );
		exit(1);
	}
	int natoms = 0;
        int step = 0;
        float time = 0;
        matrix box;
        rvec *x;
        float prec = 0;
        XDRFILE *xd;
        FILE *water, *survival, *index1, *index2 ;//, *angles ; 
	int result_xtc;  // no change till here
	int nframes=0, tot_frames=0, save_freq;
        int nh2o=17448, t0;
        int start, i, inox=0, iprot=0, axis, broke=0, dt, maxdt=300;
        int lines=0, num_water=0,num=0, indices, j; 
	float upper_cut=0.0, distsq=0.0, del=0.0 ;
        char c;
	int **resid ;
//	int C_ind[8]={0, 4, 5, 8, 10, 14, 18, 22};

        printf("THIS PROGRAM NEEDS trajectory of protein-residues without H and OW only\n");
        printf("RESIDUE-wise survival time calculation of water \n") ;

        printf("Enter resid for the survival time caln??\n");
	if (scanf("%d", &lines)==1){
        printf("Resid entered is= %d\n",lines); }
    	else {
        printf("Failed to read integer.\n");
    	}

        printf("nh2o and total frames=??\n");
        if (scanf("%d", &nh2o)==1){
        printf("nh2o=%d\n", nh2o); }
    	else {
        printf("Failed to read integer.\n");
    	}

	if (scanf("%d", &tot_frames)==1){
        printf("%d\n", tot_frames); }

	static int wat[50002][7000];
	int count[maxdt], ind[166] ;

        printf("nstart of h2os=??\n");
	if (scanf("%d", &start)==1){
        printf("%d\n",start); }

	for(i=0; i<maxdt; ++i){
	count[i]=0 ;
	}

	for(i=0; i<166; ++i){
		ind[i]=0 ;
	}
	for(t0=0; t0<tot_frames; t0++){
     		for (i=0; i<nh2o; i++) wat[t0][i]=0 ;
        }

        printf("Enter cutoff\n"  );
        if (scanf("%f", &upper_cut)==1){
        printf("cutoff=%f\n", upper_cut);}
	else {
        printf("Failed to read cutoff\n");
        }

	index1 = fopen("resid-wat-no_H_2","r");
        index2 = fopen("resid-wat-no_H_2","r");

        resid=malloc(lines * sizeof(int*));
        for (int i = 0; i < lines; ++i) {
                num=0;
        for (c = getc(index1); c != '\n'; c = getc(index1)) {
        if(fscanf(index1, "%d", &indices)==1)
                num++ ; }
		ind[i]=num ;
            resid[i]=(int*)malloc(num * sizeof(int));
        }

        for (i=0; i<lines; ++i){
	  // printf("%d\n", ind[i]+1);
           for (j=0; j<ind[i]; ++j){
             if(fscanf(index2, "%d", &indices)==1) {
               resid[i][j]=indices-1 ; // C indexing start from 0
	       //printf("resid[%d][%d]=%d\n", i, j, resid[i][j] );
             }
         else {
                 printf("Failed to read integer.\n");
             }
           }
        }

        survival = fopen("survival-time-water-around_residue","w");
        water = fopen("water-around_residue-","w");
	fprintf(survival, "#Within = %f \n", upper_cut) ;
        fprintf(water, "#Within = %f resid=%d\n", upper_cut, lines) ;
	upper_cut=upper_cut*upper_cut ;
        printf("cutoff square =%f\n", upper_cut);

        printf("Enter saving frequency of trajectory\n"  );
        if (scanf("%d", &save_freq)==1){
        printf("save_freq=%d\n", save_freq);}
	else {
        printf("Failed to read cutoff\n");
        }
	
	printf("Calling xdrfile library ...\n");
	// Opening xdrfile (xtc/trr file) with filename argv[1]
        xd = xdrfile_open(argv[1],"r");
		
        if (NULL == xd)
        	printf("ERROR: Cannot open xdrfile for reading!\n");
	// Using read_xtc_natoms function to read number of atoms (natoms)
        result_xtc = read_xtc_natoms(argv[1], &natoms);
        if (exdrOK != result_xtc)
  	      printf("ERROR: In read_xtc_natoms() \n");
        printf("Number of atoms = %d\n", natoms);
	/* Allocating N (natoms) spaces of size x[0] (DIM=3)
	x is a rvec* for 2D array of 3*natoms dimension */
	x = (rvec *)calloc(natoms, sizeof(x[0]));
	// Start looping through frames until EOF
	while(1) {
		result_xtc = read_xtc(xd, natoms, &step, &time, box, x, &prec);

		if (result_xtc == 0) { //if not EOF, read
          		if (exdrOK != result_xtc)
  	      			printf("ERROR: In read_xtc() \n");
                // code starts from here
		num_water=0;
		for(inox=start; inox<(nh2o+start); inox++)
		{ 
                   i=inox-start;
		for(iprot=0; iprot<ind[lines-1]; iprot++)
                {  
		    //Distance calculation
                      distsq = 0.0;
                      for(axis = 0; axis <= 2; axis++) {
                        del = x[inox][axis] - x[resid[lines-1][iprot]][axis];
                        if(abs(del) > box[axis][axis]*0.5){
                           del = del - copysign(box[axis][axis],del); }
                        distsq += del * del;
                      } // ends axis=0,2
                  if( distsq <= upper_cut) {
			wat[nframes][i]=1 ;
			num_water++;
			break; // breaks iprot loop and takes the next inox value
                  }
		 }//ends iprot
		} //ends inox loop
		//nframes=nframes+1 ;
		fprintf(water, "%d %d\n", nframes, num_water);
//		printf("%d %d\n", nframes, num_water);
		nframes=nframes+1 ;
		} // if (result_xtc == 0)
          	else // EOF, exit infinite loop
            		break;
	}
	// Making survival time histogram
	
//	printf("nframes=%d\n", nframes);
        for(i=0; i<nh2o; ++i){
          broke = 0;
          for(t0=0; t0 < nframes; t0++) {
              if (broke == 1) break ; // if 1; then goes to next i and breaks t0 loop
              for(dt=0; dt<=maxdt; ++dt) {
                if((t0+dt) > nframes){broke=1; break ; }

                if((wat[t0][i] == 1) && (wat[t0+dt][i] == 1))
                {
                  count[dt]++ ; // divide by bound water molecule
                 }
                else {
                        break ; // next t0 value // break dt loop
                }
             } //ends dt = 0,maxdt
          } // t0=0,nframes
        } // ends i = o ,nh2o

        // Averaging
        for(dt=0; dt<maxdt ; ++dt) {
        if(count[dt] != 0) fprintf(survival, "%d\t%f\n", dt*save_freq, (float)count[dt]/(float)(count[0])) ;
       }
 return 0;
}
