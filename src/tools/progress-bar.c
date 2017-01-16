#include "progress-bar.h"
void displayProgress(int progress){
  // int res = setupterm(NULL,fileno(stdout),NULL);
  // int row_count = tigetnum("lines");
  // int col_count = tigetnum("cols");//For terminal lines


  int col_count = 65;
  int k = 0;
  printf("\r");//将当前行全部清空，用以显示最新的进度条状态
  int j = 0;
  printf("[");

  int cols_finish = int(double(progress)/100.0*(col_count-7));

  for (j=0; j<cols_finish; j++)
          printf("=");//打印进度条上已经完成的部分，用‘+’表示
  for (j=1; j<=(col_count-7-cols_finish); j++)
          printf("_");//打印进度条上还有多少没有完成的
  fprintf(stdout, " ]%3d%%",progress);
  fflush(stdout);
  usleep(10000);
  if(progress==100)printf("\n");
}
