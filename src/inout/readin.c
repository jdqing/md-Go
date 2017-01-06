#include "readin.h"

SConfig *pSimuCfg;

bool get_words_after_keywords(char *get_words, const char *keywords, FILE *fp){
  long position, finish;
  char str[100];
  fseek( fp, 0L, SEEK_END);
  finish = ftell(fp);
  fseek( fp, 0L, SEEK_SET);
  do{
    fscanf(fp, "%s", str);
    if(strcmp(str, keywords)==0){
      fscanf(fp, "%s", get_words);
      return 1;
    }
    position = ftell(fp);
  }while( position < finish);

  return 0;
}
int get_long(FILE *fp, const char *flag, unsigned long *num){
  char str[100];
  if(get_words_after_keywords(str, flag, fp)){
    *num = (unsigned long)strtoul( str , NULL, 0);
    return 1;
  }
  return 0;
}
int get_dou (FILE *fp, const char *flag, double        *num){
  char str[100];
  if(get_words_after_keywords(str, flag, fp)){
    *num = atof( str );
    return 1;
  }
  return 0;
}
int get_str (FILE *fp, const char *flag, char       *strget){
  char str[100];
  if(get_words_after_keywords(str, flag, fp)){
    strcpy(strget,str);
    return 1;
  }
  return 0;
}
int get_bool(FILE *fp, const char *flag, bool           *bl){
  char str[100];
  if(get_words_after_keywords(str, flag, fp)){
    if(strcmp(str,"YES")==0||strcmp(str,"yes")==0||strcmp(str,"Yes")==0||
    strcmp(str,"Y")==0||strcmp(str,"y")==0){
      *bl = 1;
      return 1;
    }
    else if(strcmp(str,"NO")==0||strcmp(str,"no")==0||strcmp(str,"N")==0||
    strcmp(str,"n")==0){
      *bl = 0;
      return 1;
    }
  }
  return 0;
}
int read_Simu_Config(FILE *fp, SConfig *pCfg)
{
  get_str (fp, "PDBID",        pCfg->PDBID);
  get_dou (fp, "TemperatureK",&pCfg->TemperatureK);
  get_dou (fp, "LJScaleN",    &pCfg->LJScaleN);
  get_dou (fp, "CUTOFF",      &pCfg->CUTOFF);
  get_long(fp, "SEED",        &pCfg->SEED);

  // printf("/////////////%s\n", pCfg->PDBID);
  // printf("////////////%lf\n", pCfg->TemperatureK);
  // printf("////////////%lf\n", pCfg->LJScaleN);
  // printf("////////////%lf\n", pCfg->CUTOFF);
  // printf("////////////%ld\n", pCfg->SEED);

  return 1;
}
