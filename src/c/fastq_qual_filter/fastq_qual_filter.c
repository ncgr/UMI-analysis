/* *
 * Filters a fastq file based on the average quality of a substring
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <bstrlib.h>


// * The [libbstring] (https://github.com/msteinert/bstring) shared object library needs 
// to be in your LD\_LIBRARY\_PATH. This library can speed up your string manipulations.
// Try: LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sw/libs/bstring/0.1.1/lib/


// Could go in a header file
struct fastq_block 
{
  int header1_line;
  int sequence_line;
  int header2_line;
  int qualities_line;
  bstring header1;
  bstring sequence;
  bstring header2;
  bstring qualities;
};

// Could go in a header file
void print_block (struct fastq_block *b);
int quality (struct fastq_block *b, int offset, int length);

int main (int argc, char *argv[]) 
{
  if (argc != 8)
  {
     fprintf (stderr, "\nUsage: fastq_qual_filter infile goodfile badfile logfile min_ok offset length\n\n");
     fprintf (stderr, "   infile: path to a fastq file\n");
     fprintf (stderr, "   goodfile: path to output fastq file for passing reads\n");
     fprintf (stderr, "   badfile: path to output fastq file for failing reads\n");
     fprintf (stderr, "   logfile: path to file for logging (not used)\n");
     fprintf (stderr, "   min_ok: pass/fail mean quality threshold\n");
     fprintf (stderr, "   offset: 0-based position to begin the substring to examine\n");
     fprintf (stderr, "   length: length of substring to examine\n");
     fprintf (stderr, "\nWarning: not much arg checking is done so unexpected behavior may be encountered. For example if you enter a substring longer than the read length, or negative values. Output files will be overwritten without warning.\n\n");
     //./fastq_qual_filter disposable.fq good.fq bad.fq log 30 0 10
     exit (-1);
  } 
  
  char *fastq_file       = argv[1];
  char *output_file_good = argv[2];
  char *output_file_bad  = argv[3];
  char *log_file         = argv[4];  
  int   min_ok           = atoi (argv[5]);
  int   offset           = atoi (argv[6]);
  int   length           = atoi (argv[7]);  

  FILE *in  = fopen (fastq_file,  "r");
  FILE *out_good = fopen (output_file_good, "w");
  FILE *out_bad = fopen (output_file_bad, "w");
  FILE *log = fopen (log_file,    "w");
  
  struct fastq_block block;
  bstring line;
  block.header1_line = 0;
  block.sequence_line = 0;
  block.header2_line = 0;
  block.qualities_line = 0;  
  int block_line_counter = 0;
  int total_counter = 0; 
  int good_count = 0;
  int bad_count = 0;
  
  if (in == NULL) {return -__LINE__;}  /* some reserved var */
  
  while((line = bgets ((bNgetc) fgetc, in, (char) '\n')) != NULL)
  {
     btrimws (line);
     //printf("%s\n", line->data); 
        
     block_line_counter++;
     total_counter++;
     
     // if first line assign block.header1
     if(block_line_counter == 1)
     {
        block.header1_line = total_counter;
        block.header1 = bstrcpy(line);
        //printf("%s\n", line->data);
     }     
     // if second line assign block.sequence
     else if(block_line_counter == 2)
     {
        block.sequence_line = total_counter;
        block.sequence = bstrcpy(line);
        //printf("%s\n", line->data);
     }     
     // if third line assign block.header2
     else if(block_line_counter == 3)
     {
        block.header2 = bstrcpy(line);
        block.header2_line = total_counter;
        //printf("%s\n", line->data);
     }
     // if fourth line assign block.qualities
     else if(block_line_counter == 4)
     {
        block.qualities = bstrcpy(line);
        block.qualities_line = total_counter;
        //printf("%s\n", line->data);
               
        //print_block(&block);             
        float mean_qual = quality(&block,offset,length);
        
        if(mean_qual >= (float) min_ok)
        {
           good_count++;
           fprintf(out_good, "%s\n%s\n%s\n%s\n", block.header1->data,block.sequence->data,block.header2->data,block.qualities->data);
        }
        else
        {
           bad_count++;
           fprintf(out_bad, "%s\n%s\n%s\n%s\n", block.header1->data,block.sequence->data,block.header2->data,block.qualities->data);
        }
                     
        // reset      
        block_line_counter = 0;       
        bdestroy (block.header1);
        bdestroy (block.sequence);
        bdestroy (block.header2);
        bdestroy (block.qualities);         
     }     
     
     if( ! (total_counter % 1000000) )
     {
        fprintf (stderr, "Processed: %d lines\n", total_counter);
     }
     
     bdestroy (line);
  }
  
  fprintf (stderr, "Pass: %d\n", good_count);
  fprintf (stderr, "Fail: %d\n", bad_count);
  
  fclose(in);
  fclose(out_good);
  fclose(out_bad);
  fclose(log);
  return 0;
}


int quality (struct fastq_block *b, int offset, int length)
{
   bstring quals = bstrcpy(b->qualities);

   int total = 0;
   int i;
   
   for (i=0;  i <= length-1 ; i++)
   {
      bstring q = bmidstr (quals,0,1);
      bdelete (quals,0,1);
      char c = q->data[0];     
      int q_val = (int) c;
      q_val -= 33;      
      total += q_val;
      bdestroy(q);
   }  
   bdestroy (quals);
   return (float) (total/length);   
}

void print_block (struct fastq_block *b)
{
   //printf("%d\n", b->header1_line);     
   printf("%s\n", b->header1->data);     
   //printf("%d\n", b->sequence_line);     
   printf("%s\n", b->sequence->data);     
   //printf("%d\n", b->header2_line);     
   printf("%s\n", b->header2->data);     
   //printf("%d\n", b->qualities_line);     
   printf("%s\n\n", b->qualities->data);     
}


