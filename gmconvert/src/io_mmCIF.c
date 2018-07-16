/*

  <io_mmCIF.c>
 
  functions for reading mmCIF file into "struct BLOCK".


  LastModDate: Dec 18, 2014.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "io_mmCIF.h"




/* Functions (GLOBAL) */
void read_mmCIF_file_into_BLOCK();
void read_mmCIF_file_into_WORDs();
void translate_mmCIF_WORDs_into_BLOCK_data();
int data_malloc_from_rowline();
char* data_from_rowline();
int number_from_column();
int Integer_from_key();
struct TABLE *TABLE_from_key();
int free_WORDs();
void show_WORDs();
void free_BLOCK();
void free_ROWLINEs();
void split_into_two_words();


/* Functions (LOCAL) */
static void make_rowline_from_temp_head_word_for_one_row_table();
static struct WORD * add_new_string_to_current_end();
static void remove_quotation_from_string();
static int strip_change_line_char();
static int hash_func();
static void initialize_int_hash_tab();
static struct KEY_INT *KEY_INT_from_key();
static void insert_key_into_int_hash_tab();
static void initialize_table_hash_tab();
static struct KEY_TABLE *KEY_TABLE_from_key();
static void insert_key_into_table_hash_tab();
static void add_string_to_WORDs();
static int  number_of_WORDs();
static void free_TABLE();




void read_mmCIF_file_into_BLOCK(ifname,blk)
  char *ifname;
  struct BLOCK *blk;
{
  struct WORD HeadWord; 

  /* printf("#read_mmCIF_file_into_BLOCK('%s')\n",ifname); */
  read_mmCIF_file_into_WORDs(ifname,&HeadWord);
  /* printf("#Nword %d\n",number_of_WORDs(&HeadWord)); */
  /* show_WORDs(&HeadWord); exit(1); */
  translate_mmCIF_WORDs_into_BLOCK_data(&HeadWord,blk);
  make_rowline_from_temp_head_word_for_one_row_table(blk);
  /* printf("#block_name '%s'\n",blk->block_name); */
  /* printf("#Ntable %d\n",blk->Ntable); */
  free_WORDs(&HeadWord);
} /* end of read_mmCIF_file_into_BLOCK() */




void read_mmCIF_file_into_WORDs(ifname,HeadWord)
  char *ifname;
  struct WORD *HeadWord;
{
  FILE *fp;
  char line[MAX_LENGTH_LINE],word[MAX_LENGTH_WORD], command[MAX_LENGTH_LINE];
  char quot_state;   /* state of quotation. ';':semicolon, 's'ingle('), 'd'ouble(")  */
  int  i,L,Lline,Lword;  
  struct WORD *wn;

  /* printf("#read_mmCIF_file_into_WORDs('%s')\n",ifname); */

  L = strlen(ifname);
  if ((L>3)&&(ifname[L-3]=='.')&&(ifname[L-2]=='g')&&(ifname[L-1]=='z')){
    sprintf(command,"zcat %s",ifname);
    fp = popen(command,"r");
  }
  else{
    fp = fopen(ifname,"r");
  }
  if (fp==NULL){printf("#ERROR:Can't open ciffile \"%s\"\n",ifname); exit(1);}


  HeadWord->next = NULL;
  HeadWord->word = NULL;
  quot_state  = '-';

  wn = HeadWord;
  Lword = 0;

  while (feof(fp)==0){
    line[0] = '\0';
    fgets(line,MAX_LENGTH_LINE,fp);
    Lline = strlen(line);
    Lline = strip_change_line_char(line);
    /* printf("line '%s' quot_state '%c'\n",line,quot_state); */

    if ((Lline>0) && (line[0]!='#')){
 
      i = 0;
      while (i<Lline){

       /* semicolon (at the head of line[]) quotation */
        if ((i==0)&&(line[i]==';')&&(quot_state!='s')&&(quot_state!='d')){  
          if (quot_state==';'){quot_state = '-';}
          else                {quot_state = ';';}
        }

        /* double quotation ( the end " should be (1) at the end of line (2)the next char is space )*/
        if ((line[i]=='"')&&(quot_state!=';')&&(quot_state!='s')){ 
          if (quot_state == 'd'){
             if ((i==(Lline-1)) || (line[i+1] == ' ')) {quot_state = '-';}
          }
          else                  {quot_state = 'd';}
        }

        /* single quotation ( the end ' should be (1) at the end of line (2)the next char is space )*/
        if ((line[i]==39)&&(quot_state!='d')&&(quot_state!=';')){ 
          if (quot_state == 's'){
             if ((i==(Lline-1)) || (line[i+1] == ' ')) {quot_state = '-';}
          }
          else                  {quot_state = 's';}
        }

        /* printf("%d '%c' quot_state %c\n",i,line[i],quot_state);  */

        /** [Non-Space and Non-Quot]: extend word[] **/
        if (((line[i]!= ' ')||(quot_state==';')||(quot_state=='"')||(quot_state=='s'))){ 
          word[Lword] = line[i];
          Lword += 1;
          if (Lword>=(MAX_LENGTH_WORD-1)){
            word[Lword] = '\0';
            printf("#ERROR:Lword is over MAX_LENGTH_WORD(%d) ('%s').\n",MAX_LENGTH_WORD,word);
            exit(1);
          }
        }

        /** [Space]: terminate word[], and add to the HeadWord list. **/
        if ((line[i]== ' ')&&(quot_state=='-')){
          if (Lword>0){
            wn =  add_new_string_to_current_end(wn,word,Lword);
            Lword = 0;
          }

          while (((i+1)<Lline)&&(line[i+1]==' ')){i += 1;}
        }

        i += 1;
      }

     /** [EndOfLine and Non-Quot]: Teminate word[], and add to the HeadWord list. **/
     if ((Lword>0)&&(quot_state=='-')){
       wn = add_new_string_to_current_end(wn,word,Lword);
       Lword = 0;
     }

    }
  }

  fclose(fp);

} /* end of read_mmCIF_file_into_WORDs() */






void translate_mmCIF_WORDs_into_BLOCK_data(HeadWord,blk)
  struct WORD *HeadWord;
  struct BLOCK *blk;
{
/*
>> example of mmCIF lines << 
data_1UBQ
#
_entry.id   1UBQ               STATE 'O' 
#
_audit_conform.dict_name       mmcif_pdbx.dic                    STATE 'O'
_audit_conform.dict_version    4.007                             STATE 'O'
_audit_conform.dict_location   http://mmcif.pdb.org/dictionaries/ascii/mmcif_pdbx.dic STATE 'O'
#
_database_2.database_id     PDB            STATE 'O'
_database_2.database_code   1UBQ           STATE 'O'
#
loop_                                      STATE 'L'
_database_PDB_rev.num                      STATE 'L'
_database_PDB_rev.date                     STATE 'L'
_database_PDB_rev.date_original            STATE 'L'
_database_PDB_rev.status                   STATE 'L'
_database_PDB_rev.replaces                 STATE 'L'
_database_PDB_rev.mod_type                 STATE 'L'
1 1987-04-16 1987-01-02 ? 1UBQ 0           STATE 'R'
2 1987-07-16 ?          ? 1UBQ 1           STATE 'R'
3 2003-04-01 ?          ? 1UBQ 1           STATE 'R'
4 2009-02-24 ?          ? 1UBQ 1           STATE 'R'
5 2011-03-09 ?          ? 1UBQ 1           STATE 'R'
#
*/
  char STATE; /* 'O'ne_row, 'L'oop_item, 'R'ow_data for the loop table. */
/*
 >> rule of 'STATE' <<
  if ('loop_') -> 'L'
  if (not 'L') and ('_') -> 'O'
  if ('L') and '(not '_') -> 'R'
 >> transition rule of 'STATE' << 
  'O'->'O', 'O'->'L' 
  'L'->'L', 'L'->'R' 
  'R'->'R', 'R'->'O','R'->'L'
*/

  char category[MAX_LENGTH_LINE], item[MAX_LENGTH_LINE];
  struct WORD *wn;
  struct TABLE *tb,*tb_exist;
  struct ROWLINE *rn;
  int    c,extend;

  blk->block_name = NULL;
  blk->Ntable = 0;
  blk->head_table.next = NULL;
  initialize_table_hash_tab(blk->table_hash_tab);
  tb = &(blk->head_table);
  tb->Ncolumn = 0;

  /* printf("#translate_mmCIF_WORDs_into_BLOCK_data()\n"); */

  wn = HeadWord;
  STATE = '-';
  rn = NULL;
 
  while (wn->next != NULL){
    wn = wn->next;

/*    printf("'%s' %c\n",wn->word,STATE);  */

    /* (1) data_[PDB_ID] */
    if ((strlen(wn->word)>=5) && (strncmp(wn->word,"data_",5)==0)){
       blk->block_name = (char *)malloc(sizeof(char)*(strlen(wn->word)+1)); 
       sprintf(blk->block_name,"%s",wn->word);
    }    
    /* (2) loop_ */
    else if (strncmp(wn->word,"loop_",5)==0){
      STATE = 'L';
    }

    /* (3) '_[category].[item]' */
    else if (wn->word[0]=='_'){
      split_into_two_words(wn->word,1,'.',category,item);
      tb_exist = TABLE_from_key(blk->table_hash_tab,category);
      if (tb_exist==NULL){
        /* printf("##malloc new tb for key '%s'\n",category); fflush(stdout); */
        tb->next = (struct TABLE*)malloc(sizeof(struct TABLE));
        tb = tb->next;
        tb->Ncolumn = 0;
        tb->Nrow    = 0;
        tb->next = NULL;
        tb->table_name = (char *)malloc(sizeof(char)*(strlen(category)+1));
        sprintf(tb->table_name,"%s",category);
        tb->head_column_name.next = NULL;
        tb->head_rowline.next = NULL;
        tb->temp_head_word.next = NULL;
        rn = &(tb->head_rowline);  
        initialize_int_hash_tab(tb->colnum_hash_tab);
        insert_key_into_table_hash_tab(blk->table_hash_tab,category,tb);
        blk->Ntable += 1;
        tb_exist = tb;
      }
      add_string_to_WORDs(item,&(tb->head_column_name));
      insert_key_into_int_hash_tab(tb->colnum_hash_tab,item,tb->Ncolumn);
      tb->Ncolumn += 1;
 
      if ((wn->next != NULL) && (wn->next->word[0] != '_') && (STATE!='L')){
        STATE = 'O';
      }
    }

    /* (4) values for the loop */
    else if ((wn->word[0]!='_')&&(strlen(wn->word)>0)&&(tb->Ncolumn>0)&&((STATE =='L')||(STATE=='R'))){
       /* printf("#malloc ROWLINE with Ncolumn=%d\n",tb->Ncolumn); */
       rn->next = (struct ROWLINE*)malloc(sizeof(struct ROWLINE));
       rn = rn->next;
       tb->Nrow += 1;
       rn->next = NULL;
       rn->data = (char **)malloc(sizeof(char *)*tb->Ncolumn);
       STATE = 'R';
       c = 0;
       do{
         if (wn->word[0]=='_'){
           printf("#ERROR(line '%s' table '%s'). Not enough items(c %d) for Ncolumn=%d  in data row. \n",wn->word,tb->table_name,c,tb->Ncolumn); 
           exit(1);
         }
         remove_quotation_from_string(wn->word);
         rn->data[c] = (char *)malloc(sizeof(char)*(strlen(wn->word)+1));
         sprintf(rn->data[c],"%s",wn->word);
         /* printf(">'%s' %d '%s'\n",tb->table_name,c,rn->data[c]); */
         extend = 0;
         if ((c<(tb->Ncolumn-1))&&(wn->next != NULL)&&(wn->next->word[0]!='_')&&(strcmp(wn->next->word,"loop_")!=0)){ 
           wn = wn->next;
           extend = 1;
           c += 1;
         }
       } while ((c<tb->Ncolumn)&&(extend==1));
    }

    /* (5) values for one-row table */
    else if ((wn->word[0]!='_')&&(strlen(wn->word)>0)&&(STATE =='O')){
        /* printf("one-row values (table '%s') add_string_to_WORDs('%s');\n",tb->table_name,wn->word); */
        remove_quotation_from_string(wn->word);
        add_string_to_WORDs(wn->word,&(tb->temp_head_word));
    }
  }


} /* end of translate_mmCIF_WORDs_into_BLOCK_data() */




void make_rowline_from_temp_head_word_for_one_row_table(blk)
  struct BLOCK *blk;
{
  struct TABLE *tb;
  struct ROWLINE *rn;
  struct WORD    *ln;
  int c;
 
  tb = &(blk->head_table);
  while (tb->next != NULL){
    tb = tb->next;
    if ((tb->Nrow==0) && (tb->Ncolumn>0) && (tb->Ncolumn == number_of_WORDs(&(tb->temp_head_word)))){
      rn = &(tb->head_rowline);  
      rn->next =  (struct ROWLINE*)malloc(sizeof(struct ROWLINE));
      rn = rn->next;
      rn->next = NULL; 
      rn->data = (char **)malloc(sizeof(char *)*tb->Ncolumn);
      tb->Nrow += 1;
      ln = &(tb->temp_head_word);
      c = 0;
      while (ln->next != NULL){
        ln = ln->next;
        rn->data[c] = (char *)malloc(sizeof(char)*(strlen(ln->word)+1));
        sprintf(rn->data[c],"%s",ln->word);
        c += 1;
      }
      free_WORDs(&(tb->temp_head_word));
    }
  }
} /* end of make_rowline_from_temp_head_word_for_one_row_table() */




struct WORD * add_new_string_to_current_end(endnode,word,Lword)
  struct WORD *endnode; /* the current end node */
  char *word;
  int  Lword;
{
  int i;

  if (Lword>0){
    word[Lword] = '\0';

    /* add new struct LINE */
    endnode->next = (struct WORD*)malloc(sizeof(struct WORD));
    endnode->next->word = (char *)malloc(sizeof(char)*(Lword+1));
    /* sprintf(endnode->next->word,"%s",word); */
    for (i=0;i<=Lword;++i){
      endnode->next->word[i] = word[i];
    }
    endnode->next->word[Lword] = '\0';
    endnode->next->next = NULL;
    endnode = endnode->next;
  }
  return(endnode); /* return 'new' endnode */
}



void remove_quotation_from_string(word)
  char *word;
{
  int Lword,i;
/*
   remove the head and the tail by semicolon, double quotation, single quotation.
         for example, ;ABC; , "ABC" , 'ABC' ==> ABC, ABC, ABC.
>>example <<
 Lword = 5 

        01234
word[]= 'abc'   word[5] = '\0'
     => abc     word[3] = '\0'

*/
  Lword = strlen(word);
  if ( ((word[0]==';')&&(word[Lword-1]==';')) ||
       ((word[0]=='"')&&(word[Lword-1]=='"')) ||
       ((word[0]== 39)&&(word[Lword-1]== 39)) ){
    for (i=1;i<=(Lword-2);++i){
      word[i-1] = word[i];
    } 
    word[Lword-2] = '\0';
  }
}





int strip_change_line_char(str)
  char *str;
{
  int L;
  L = strlen(str);
  if (str[L-1]=='\n') {str[L-1] = '\0'; L -= 1;}
  if (str[L-1]=='\r') {str[L-1] = '\0'; L -= 1;}
  if (str[L-1]=='\n') {str[L-1] = '\0'; L -= 1;}
  if (str[L-1]=='\r') {str[L-1] = '\0'; L -= 1;}
  return(L);
}





void split_into_two_words(line,start_pos,splsym,first_word,second_word)
  char *line;
  int  start_pos;  /* start position of line[]. the default should be 0 */
  char splsym;     /* symbol for spliting. for example, '.', '\t',' ',.. */
  char *first_word;
  char *second_word;
{
 /*
  if line[] = '_audit_conform.dict_name'
     splsym = '.', start_pos = 1.
  ==> 
   first_word[]  = 'audit_conform'
   second_word[] = 'dict_name'
 
  if (line[] = '_audit_conform.dict_name       mmcif_pdbx.dic'
     splsym = ' ', start_pos = 1.

 ==> first_word[]  =  'audit_conform.dict_name'
     second_word[] = 'mmcif_pdbx.dic'


*/
  int i,j,L;
  char findspl;

  L = strlen(line);
  i = start_pos;
  
  findspl = 0;
  while ((i<L)&&(findspl==0)){
    if (line[i]==splsym){
      findspl = 1;
      for (j=start_pos;j<i;++j){ 
        first_word[j-start_pos] = line[j];
      } 
      first_word[i-start_pos] = '\0';
      while ((line[i]==splsym) && (i<L)){
        i += 1;
      }

      for (j=i;j<L;++j){ 
        second_word[j-i] = line[j];
      }
      second_word[L-i] = '\0';
    }
   ++i;
  }
 
  if (findspl==0){
    for (j=start_pos;j<L;++j){ 
      first_word[j-start_pos] = line[j];
    }
    first_word[L] = '\0';
    second_word[0] = '\0';
  }

} /* end of split_into_two_words() */





int hash_func(key)
  char *key;
{
  int i,s;
  s = 0;
  for (i=0;i<strlen(key);++i){
    s += key[i];
  }
  return(s % HASH_SIZE);
}

void initialize_int_hash_tab(int_hash_tab)
  struct KEY_INT *int_hash_tab;
{
  int i;
  for (i=0;i<HASH_SIZE;++i){
    int_hash_tab[i].next = NULL;
  }  
}


int data_malloc_from_rowline(new_str,tb,rn,column_name)
  char   *new_str;
  struct TABLE *tb;
  struct ROWLINE *rn;
  char *column_name;
{
  struct KEY_INT *p;
  int c,L;

  printf("#data_malloc_from_rowline(column_name '%s')\n",column_name);
  p = &(tb->colnum_hash_tab[hash_func(column_name)]);
  c = -1;
  while (p->next !=NULL){
    p = p->next;
    if (strcmp(column_name,p->key)==0){
      c = p->int_value; 
    }
  } 
  
  if ((c<0)||(c>tb->Ncolumn)){ return(0); }
  
  else{
    L = strlen(rn->data[c]);
    new_str = (char *)malloc(sizeof(char)*(L+1));
    sprintf(new_str,"%s",rn->data[c]);
    printf("#data_malloc_from_rowline(column_name '%s' new_str '%s')\n",column_name,new_str);
    return(L);
  }
}


char* data_from_rowline(tb,rn,column_name)
  struct TABLE *tb;
  struct ROWLINE *rn;
  char *column_name;
{
  struct KEY_INT *p;
  int c;

  p = &(tb->colnum_hash_tab[hash_func(column_name)]);
  c = -1;
  while (p->next !=NULL){
    p = p->next;
    if (strcmp(column_name,p->key)==0){
      c = p->int_value; 
    }
  } 
  
  if ((c<0)||(c>tb->Ncolumn)){ return(""); }
  
  else{
    return(rn->data[c]);
  }
}




int number_from_column(tb,column_name)
  struct TABLE *tb;
  char *column_name;
{
  struct KEY_INT *p;
  p = &(tb->colnum_hash_tab[hash_func(column_name)]);
  while (p->next !=NULL){
    p = p->next;
    if (strcmp(column_name,p->key)==0){
      return(p->int_value); 
    }
  } 
  return(-1);
} /* end of number_from_column() */



int Integer_from_key(hash_table,key)
  struct KEY_INT *hash_table;
  char *key;
{
  struct KEY_INT *p;

  p = &(hash_table[hash_func(key)]);
  while (p->next !=NULL){
    p = p->next;
    /* printf("query %s qhash %d key '%s' value '%s'\n",key,hash(key),p->key,p->value); */
    if (strcmp(key,p->key)==0){
      return(p->int_value); 
    }
  } 

  return(-1);

} /* end of Integer_from_key() */






struct KEY_INT *KEY_INT_from_key(hash_table,key)
  struct KEY_INT *hash_table;
  char *key;
{
  struct KEY_INT *p;

  p = &(hash_table[hash_func(key)]);
  while (p->next !=NULL){
    p = p->next;
    if (strcmp(key,p->key)==0){ return(p); }
  } 
  return(NULL);
}



void insert_key_into_int_hash_tab(int_hash_tab,key,value)
  struct KEY_INT *int_hash_tab;
  char *key;
  int  value;
{
  struct KEY_INT *p,*last_cell;
  int h;

/*
  printf("#insert_key_into_int_hash_tab(key '%s' hash %d value %d)\n",key,hash_func(key),value);
*/  
  p = KEY_INT_from_key(int_hash_tab,key);

  if (p != NULL){
    /* If KEY_INT with key exists, overwrite existing KEY_INT. */
    p->int_value = value;
  }
  else{
    /* If KEY_INT with key does not exist,  malloc new KEY_INT  */
    p = (struct KEY_INT *)malloc(sizeof(struct KEY_INT));

    if (p==NULL){
      printf("#ERROR:out of memory.\n");
      exit(1);
    }

    p->key  = (char *)malloc(sizeof(char)*(strlen(key)+1));
    sprintf(p->key, "%s",key);
    p->int_value = value;
    p->next = NULL;

    h = hash_func(key);
    last_cell = int_hash_tab[h].next;
    p->next = last_cell;
    int_hash_tab[h].next = p;
  }
}





void initialize_table_hash_tab(table_hash_tab)
  struct KEY_TABLE *table_hash_tab;
{
  int i;
  for (i=0;i<HASH_SIZE;++i){
    table_hash_tab[i].next = NULL;
  }  
}


struct TABLE *TABLE_from_key(hash_table,key)
  struct KEY_TABLE *hash_table;
  char *key;
{
  struct KEY_TABLE *p;

/*
  printf("#TABLE_from_key(key '%s' hash_func %d)\n",key,hash_func(key));
*/

  p = &(hash_table[hash_func(key)]);
  while (p->next !=NULL){
    p = p->next;
    if (strcmp(key,p->key)==0){ return(p->table_value); }
  } 
  return(NULL);
}



struct KEY_TABLE *KEY_TABLE_from_key(hash_table,key)
  struct KEY_TABLE *hash_table;
  char *key;
{
  struct KEY_TABLE *p;

  p = &(hash_table[hash_func(key)]);
  while (p->next !=NULL){
    p = p->next;
    if (strcmp(key,p->key)==0){ return(p); }
  } 
  return(NULL);
}





void insert_key_into_table_hash_tab(table_hash_tab,key,table)
  struct KEY_TABLE *table_hash_tab;
  char *key;
  struct TABLE *table;
{
  struct KEY_TABLE *p,*last_cell;
  int h;
/*
  printf("#insert_key_into_table_hash_tab(key '%s' hash %d)\n",key,hash_func(key));
 */ 
  p = KEY_TABLE_from_key(table_hash_tab,key);

  if (p != NULL){
    /* If KEY_INT with key exists, overwrite existing KEY_INT. */
    p->table_value = table;
  }
  else{
    /* If KEY_INT with key does not exist,  malloc new KEY_INT  */
    p = (struct KEY_TABLE *)malloc(sizeof(struct KEY_TABLE));

    if (p==NULL){
      printf("#ERROR:out of memory.\n");
      exit(1);
    }

    p->key  = (char *)malloc(sizeof(char)*(strlen(key)+1));
    sprintf(p->key, "%s",key);
    p->table_value = table;
    p->next = NULL;

    h = hash_func(key);
    last_cell = table_hash_tab[h].next;
    p->next = last_cell;
    table_hash_tab[h].next = p;
  }
}







void add_string_to_WORDs(newstr,HeadWordNode)
 char *newstr;
 struct WORD *HeadWordNode;
{
  struct WORD *newN,*lastN;

  /* printf("# add_string_to_WORDs(newstr '%s')\n",newstr); */
  lastN = HeadWordNode;

  while (lastN->next != NULL){
    lastN = lastN->next;
  }

  lastN->next = (struct WORD *)malloc(sizeof(struct WORD));
  newN = lastN->next;
  newN->next = NULL;
  newN->word = (char *)malloc(sizeof(char)*(strlen(newstr)+1));
  sprintf(newN->word,"%s",newstr);

} /* end of add_string_to_WORDs() */




/*
int free_WORDs_orig(HeadWordNode)
 struct WORD *HeadWordNode;
{
 struct WORD *ln;
 int nfree;

 nfree = 0; 

 ln = HeadWordNode;
 while (ln->next != NULL){ln = ln->next;}

 while ((ln->prev!=NULL)&&(ln->prev!=HeadWordNode)){
   ln = ln->prev;
   if (ln->next != NULL){
     if (ln->next->word != NULL) {free(ln->next->word);}
     free(ln->next);
     nfree += 1;
   }
 }

 HeadWordNode->next = NULL;
 HeadWordNode->word = NULL;
 printf("#nfree %d\n",nfree);
 return(nfree);
}
*/




int free_WORDs(HeadWordNode)
 struct WORD *HeadWordNode;
{
 struct WORD *current, *next;
 int nfree;

 nfree = 0; 
 current = HeadWordNode;
 next = current->next;
 while (next != NULL){
   current = next;
   next = current->next;
   free(current->word);
   current->word = NULL;
   free(current);
   current = NULL;
   nfree += 1; 
 } 

 HeadWordNode->next = NULL;
 HeadWordNode->word = NULL;
/*
 printf("#nfree %d\n",nfree);
*/
 return(nfree);

} /* end of free_WORDs() */



void show_WORDs(HeadWordNode)
 struct WORD *HeadWordNode;
{
  struct WORD *bn;
  bn = HeadWordNode;
  while (bn->next != NULL){
    bn = bn->next;
    printf("'%s'\n",bn->word);
  }
} /* end of show_WORDs() */


int  number_of_WORDs(HeadWordNode)
 struct WORD *HeadWordNode;
{
  struct WORD *bn;
  int n;
  n = 0;
  bn = HeadWordNode;
  while (bn->next != NULL){
    bn = bn->next;
    n += 1;
  }
  return(n);
} /* end of number_of_WORDs() */



void free_ROWLINEs(head_rowline,Ncolumn)
  struct ROWLINE *head_rowline;
{
 struct ROWLINE *current, *next;
 int i;

 current = head_rowline;
 next = current->next;
 while (next != NULL){
   current = next;
   next = current->next;
   for (i=0;i<Ncolumn;++i){
     free(current->data[i]);
   } 
   free(current->data);
   current->data = NULL;
   free(current);
 } 

 head_rowline->next = NULL;
 head_rowline->data = NULL;

} /* end of free_ROWLINEs() */


void free_BLOCK(blk)
  struct BLOCK *blk;
{
  struct TABLE     *tb_curr,*tb_next;
  struct KEY_TABLE *kt_curr,*kt_next;
  int i;
 

  tb_curr = &(blk->head_table);  
  tb_next = tb_curr->next;
  while (tb_next != NULL){
    tb_curr = tb_next;
    tb_next = tb_curr->next;
    free_TABLE(tb_curr);
  } 

  for (i=0;i<HASH_SIZE;++i){
    kt_curr = &(blk->table_hash_tab[i]);
    kt_next = kt_curr->next;
    while (kt_next != NULL){
      kt_curr = kt_next;
      kt_next = kt_curr->next;
      free(kt_curr);
    } 
  }

  free(blk->block_name);

} /* end of free_BLOCK() */


void free_TABLE(tb)
  struct TABLE *tb;
{
  int i;
  struct KEY_INT *ki_curr,*ki_next;

/*
  printf("#free_TABLE('%s')\n",tb->table_name);
*/


  free_ROWLINEs(&(tb->head_rowline),tb->Ncolumn);
  free_WORDs(&(tb->head_column_name));

  for (i=0;i<HASH_SIZE;++i){
    ki_curr = &(tb->colnum_hash_tab[i]);
    ki_next = ki_curr->next;
    while (ki_next != NULL){
      ki_curr = ki_next;
      ki_next = ki_curr->next;
      free(ki_curr);
    } 
  }

  free_WORDs(&(tb->temp_head_word));
  free(tb->table_name);

  tb->Ncolumn = 0;
  tb->Nrow    = 0;
} /* end of free_TABLE() */

