/* Operations on bits, 
To be used by decode_mama

JLP
Version 08-10-91

main()
{
short int i1, i2, i3;

i1 = 1; i2 = 3; 
jlp_orshort(&i1,&i2,&i3);
printf(" %d or %d is %d :\n",i1,i2,i3);

i1 = 4; i2 = 7; 
jlp_andshort(&i1,&i2,&i3);
printf(" %d and %d is %d :\n",i1,i2,i3);

i1 = 4; i2 = 2; 
jlp_shftshort(&i1,&i2,&i3);
printf(" %d shifted of %d is %d :\n",i1,i2,i3);

i1 = 4; i2 = 2; 
jlp_btestshort(&i1,&i2,&i3);
printf(" bit # %d of %d is %d :\n",i2,i1,i3);

i1 = 4; i2 = 4; 
jlp_btestshort(&i1,&i2,&i3);
printf(" bit # %d of %d is %d :\n",i2,i1,i3);

}
#define DEBUG 1
*/

int jlp_orshort(i1,i2,i3)
short int *i1, *i2, *i3;
{
 *i3 = *i1 | *i2;
#if DEBUG
 printf(" %d or %d is %d :\n",*i1,*i2,*i3);
#endif
 return (0);
}
int jlp_andshort(i1,i2,i3)
short int *i1, *i2, *i3;
{
 *i3 = *i1 & *i2;
#if DEBUG
 printf(" %d and %d is %d :\n",*i1,*i2,*i3);
#endif
 return (0);
}
int jlp_shftshort(i1,i2,i3)
short int *i1, *i2, *i3;
{
 short int i;

 if(*i2 > 0)
   *i3 = *i1 << *i2;
 else
   {i = -1 * *i2;
   *i3 = *i1 >> i;}

#if DEBUG
 printf(" %d shifted of %d is %d :\n",*i1,*i2,*i3);
#endif
 return (0);
}
int jlp_btestshort(i1,i2,i3)
short int *i1, *i2, *i3;
{
 short int itest;

 itest = *i1 & ( 1 << *i2);
 if(itest == 0) *i3 = 0;
 else *i3 = 1;

#if DEBUG
 printf(" bit number %d of %d is %d :\n",*i2,*i1,*i3);
#endif

 return (0);
}
