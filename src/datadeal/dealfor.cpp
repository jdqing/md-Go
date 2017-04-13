#include "dealfor.h"
const int all=21000000;
double R[all];
double radius_of_gyration ( int n1, int n2 )
{
	int i;

	double xx [3], r1 = 0;

	double xcm = 0, ycm = 0, zcm = 0;

	int N = 0;

	////////////////////////////////////////////////////

	for (i = n1; i < n2 + 1; i++) { N++; xcm += x [i]; ycm += y [i]; zcm += z [i]; }

	xcm /= N; ycm /= N; zcm /= N;

	////////////////////////////////////////////////////

	for (i = n1; i < n2 + 1; i++)
	{
		xx [0] = x [i] - xcm;

		xx [1] = y [i] - ycm;

		xx [2] = z [i] - zcm;

		r1 += xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];
	}

	r1 = sqrt ( r1 / N );

	return ( r1 );
}

double remaining_native_contacts ()
{
//	FILE * f1;

	int i, j, k, j1;

	double count = 0;

	double xx [3], r2;

	///////////////////////////////////////////////////////

	for (k = 1; k < native_list_mass + 1; k += 2)
	{
		j1 = (k + 1) / 2;

		///////////////////////////////////////////////////////

		i = native_list_content [k];

		j = native_list_content [k + 1];

		///////////////////////////////////////////////////////

		xx [0] = x [i] - x [j];

		xx [1] = y [i] - y [j];

		xx [2] = z [i] - z [j];

		r2 = xx [0] * xx [0] + xx [1] * xx [1] + xx [2] * xx [2];

		///////////////////////////////////////////////////////

		if ( r2 < D2_LJ_n [j1] )	count += 1.0;

		else if ( r2 < D2_LJ_CUTOFF_n [j1] )
		{
			xx [0] = D2_LJ_n [j1] / r2;

			xx [1] = xx [0] * xx [0];

			xx [2] = xx [1] * xx [0];

			xx [0] = xx [2] * xx [2];

			count += ( 2.0 * xx [2] - xx [0] );
		}


		else count += 0;
	}

	///////////////////////////////////////////////////////

	return ( count / N0_NATIVE );
}
void Caculate_TA_MSD()
{

	int delta[]={1,2,3,4,5,6,7,8,9,
							 1e1,2e1,3e1,4e1,5e1,6e1,7e1,8e1,9e1,
							 1e2,2e2,3e2,4e2,5e2,6e2,7e2,8e2,9e2,
						   1e3,2e3,3e3,4e3,5e3,6e3,7e3,8e3,9e3,
					     1e4,2e4,3e4,4e4,5e4,6e4,7e4,8e4,9e4,
						   1e5,2e5,3e5,4e5,5e5,6e5,7e5,8e5,9e5};
	int t[]    ={2e1,2e3,1e5,2e6};
	//
	int lengthdelta = sizeof(delta)/sizeof(delta[0]);
	int lengtht     = sizeof(t)/sizeof(t[0]);

	long kkk;
	int i,j,k,l,m,countt;
	//
	double qq;

	char outfile1[]="NAF332.dat";
	char outfile2[]="TA-MSD332.dat";
	char openfile[]="Coordinates_1PIN_T332.dat";
	FILE *fp;
	fp = fopen (openfile,"r");
	if (fp!=NULL)
	{
		double dxx,dyy,dzz;
		long position,finish;
		fseek(fp,0L,SEEK_END);
		finish=ftell(fp);
		fseek(fp,0L,SEEK_SET);
		countt=0;
		for (j = 0; j < all; j++) {
			// position=ftell(fp);
			// if(position>=finish)break;

			fread(x,sizeof(double),NTP_atom+1,fp);
			fread(y,sizeof(double),NTP_atom+1,fp);
			fread(z,sizeof(double),NTP_atom+1,fp);
			dxx=x[NTP_atom]-x[1];
			dyy=y[NTP_atom]-y[1];
			dzz=z[NTP_atom]-z[1];
			R[j]=sqrt(dxx*dxx+dyy*dyy+dzz*dzz);
			// fscanf(fp,"%ld %le %le\n", &kkk,R+j,&qq);
			if(j%100000==0)
			printf("%d	%lf\n", j,R[j] );
		}
		countt=j;
		printf("%d\n",countt );

		fclose (fp);
	}
	// char openfile[]="Dimensions_1PIN_T332.dat";
	// FILE *fp;
	// fp = fopen (openfile,"r");
	// if (fp!=NULL)
	// {
	// 	long position,finish;
	// 	fseek(fp,0L,SEEK_END);
	// 	finish=ftell(fp);
	// 	fseek(fp,0L,SEEK_SET);
	// 	countt=0;
	// 	for (j = 0; j < all; j++) {
	// 		position=ftell(fp);
	// 		if(position>=finish)break;
	// 		fscanf(fp,"%ld %le %le\n", &kkk,R+j,&qq);
	// 		// printf("%ld	%lf	%lf\n", kkk,R[j],qq );
	// 	}
	// 	countt=j;
	// 	printf("%d\n",countt );
	//
	// 	fclose (fp);
	// }


	double sum1,av1,sum2,av2,sum3,av3,sum4,av4;

	fp = fopen (outfile1,"w");
	fclose(fp);
	fp = fopen (outfile1,"a+");

	for(i=0;i<lengthdelta;i++)
	{
		fprintf(fp,"%d",delta[i] );
		for(l=lengtht-1;l>=0;l--)
		{
			if(t[l]<2*delta[i])continue;
			if(t[l]>countt)continue;
			int avnum;
			avnum = countt/t[l];
			// printf("%d\n",avnum );
			sum1=0;
			for(m=0;m<avnum;m++)
			{
				sum3=0;
				for(j=m*t[l];j<(m+1)*t[l];j++)
				{
					sum3 = sum3+R[j];
				}
				av3=sum3/t[l];

				sum2 = 0;
				for(j=m*t[l];j<(m+1)*t[l]-delta[i];j++)
				{
					// sum2 =sum2 + (R[j+delta[i]]-R[j])*(R[j+delta[i]]-R[j]);TA-MSD
					sum2 =sum2 + (R[j]-av3)*(R[j+delta[i]]-av3);
				}
				av2=sum2/(t[l]-delta[i]);

				sum4=0;
				for(j=m*t[l];j<(m+1)*t[l];j++)
				{
					sum4=sum4+(R[j]-av3)*(R[j]-av3);
				}
				av4=sum4/t[l];
				double tmp;
				tmp=av2/av4;

				sum1=sum1+tmp;
			}
			av1=sum1/avnum;
			if(av1<=0)continue;
			fprintf(fp,"		%lf",av1 );
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	fp = fopen (outfile2,"w");
	fclose(fp);
	fp = fopen (outfile2,"a+");


	for(i=0;i<lengthdelta;i++)
	{
		fprintf(fp,"%d",delta[i] );
		for(l=lengtht-1;l>=0;l--)
		{
			if(t[l]<2*delta[i])continue;
			if(t[l]>countt)continue;
			int avnum;
			avnum = countt/t[l];
			// printf("%d\n",avnum );
			sum1=0;
			for(m=0;m<avnum;m++)
			{
				sum2 = 0;
				for(j=m*t[l];j<(m+1)*t[l]-delta[i];j++)
				{
					sum2 =sum2 + (R[j+delta[i]]-R[j])*(R[j+delta[i]]-R[j]);
				}
				av2=sum2/(t[l]-delta[i]);

				sum1=sum1+av2;
			}
			av1=sum1/avnum;
			fprintf(fp,"		%lf",av1 );
		}
		fprintf(fp,"\n");
	}
	fclose(fp);


}

int deal_for( const char *flag)
{


	// Caculate_TA_MSD_ZWT();

	return 0;

}
