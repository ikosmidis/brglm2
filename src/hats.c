void hatsc (int* n, int* p, double* x, double* invfisherinf, double* w, double* hat)
{
  int i, j, o, pos1, pos2, pos3;
  double summi ; /*summ[*n],*/
  for (o=0;o<*n;o++) {
    /*summ[o] = 0 ;*/
    summi = 0;
    pos1 = o * *p ;
    for (i=0;i<*p;i++) {
      pos2 = i * *p;
      pos3 = pos1 + i;
      for (j=0;j<*p;j++) {
	/*summ[o] += x[o * *p + i]*x[o * *p + j]*invfisherinf[j + *p *i]; */
	/*summi += x[o * *p + i]*x[o * *p + j]*invfisherinf[j + *p *i];*/
	summi += x[pos3]*x[pos1 + j]*invfisherinf[j + pos2];
      }
    }
    /*hat[o] = summ[o] * w[o];*/
    hat[o] = summi * w[o];
  }
}
