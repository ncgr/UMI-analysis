/*
 * Feed on a list of positions, UMI counts, and their UMIs
 * Output UMIs and all of their locations
 * Must use output of umi_counterizer.pl
 *
 *     Chr01:21772 6 AAGGGCTCTG TCTCTCTGGC CGTTCACAGG TGATTCGATG GGTTTTTAGG AGTTCAGTCG
 * etc.
 */
  

/*
 * Using UTHASH to record which UMIs have been used to find all of their positions. It is efficient and
 * fast.  
 * 
 * http://troydhanson.github.io/uthash/userguide.html
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "/path/to/uthash.h"

//
// Don't forget to edit these to match your umi counts file
//

#define UMI_WIDTH        22
#define POSITION_STR_LEN 64
#define LINE_BYTES       1000000

/*
 * This is the prototype of a hashable struct that maps UMIs to counts
 */
struct count_strings 
{
   char umi[UMI_WIDTH];     // The key
   int count;               // The value
   UT_hash_handle hh;       // This makes it hashable
};

/*
 * Prototype of the string sorting function
 */
int key_sort(struct count_strings *a, struct count_strings *b);

int main (int argc, char *argv[]) 
{
   if (argc < 2)
   {
      fprintf (stderr, "\n Usage: umi_analyzer inputfile ...\n\n") ;
      fprintf (stderr, "\n e.g:  ./exe /home/projects/adhoc/vulcanite/run_170420/star_output_intronless/umi_counts_5000.txt...\n\n") ;
      exit (-1);
   }
 
   char *file          = argv[1];
   FILE *fi = fopen(file, "rb");
   FILE *fj = fopen(file, "rb");
   
   if (fj == NULL) 
   {
        fprintf(stderr, "Unable to open %s", file);
        exit(-1);
   }
   
   if (fi == NULL) 
   {
        fprintf(stderr, "Unable to open %s", file);
        exit(-1);
   }

   // This is the head of the hash of used UMIs
   struct count_strings *umi_counts = NULL; 

   // Used by strtok_r in i and j loops
   const char sep[2]   = " ";  
   
   // Outer loop vars
   int  i_line_counter  = 1;
   
   // Get number of lines in the input file
   int num_lines = get_num_lines(file);
   
   char i_line_str[LINE_BYTES];
   int  i_bytes_to_read = LINE_BYTES; // Arbitrary... Needs to be longer than the longest line in the file.
   while(fgets (i_line_str, i_bytes_to_read, fi) != NULL)
   { 
      if(!(i_line_counter % 1000))
      {
          fprintf(stderr,"Processed %d of %d lines \n",i_line_counter,num_lines);
      }

      char i_position[POSITION_STR_LEN];
      int  i_umi_count     = 0;
      int  i_token_counter = 1;
      char * itoken        = NULL;   
      int i_byte = ftell(fi);
      
      // Eliminate newline
      char *newline = strchr( i_line_str, '\n' );
      if ( newline )
      *newline = 0;
      
      // This gets tokens of the outer string
      char * iptr;
      itoken = strtok_r(i_line_str, sep, &iptr);
           
      while( itoken != NULL )
      {        
         if(i_token_counter == 1)
         {
            strcpy(i_position,itoken);
         }
         else if(i_token_counter == 2)
         {
            i_umi_count = atoi(itoken);
         }
         else
         {
            // i_token_counter must be > 2 and itoken is therefore a UMI
            //printf( "I: %d   %s   %d   %s\n", i_token_counter,i_position,i_umi_count,itoken );
            
            // Check whether this UMI has been searched with before
            struct count_strings *umi_searched_already;

            // See if this UMI is NOT in the hash
            HASH_FIND_STR(umi_counts,itoken,umi_searched_already);
                         
            //UMI used before             
            if(umi_searched_already != NULL) 
            {
               //printf( "Found: %s\n", itoken );
               
               //
               // The flow must go to the next iteration of the while itoken
               //
            }
            // New UMI
            else                            
            {
               //printf( "New:   %s\n", itoken );
               
               // Add  UMI with count of 1 to the hash
               struct count_strings *s;
               s = malloc(sizeof(struct count_strings));
               strcpy(s->umi,itoken);
               s->count = 1;
               HASH_ADD_STR(umi_counts,umi,s);
               
               // We are at the current i UMI (itoken). We want a list of positions for that UMI
               // starting with the outer string i position
               //
               // MAIN OUTPUT
               //
               printf("%s %s ",itoken,i_position);
                                            
               int j_line_counter = 1;
               
               // Go to the line after the current i_line
               fseek(fj,i_byte,SEEK_SET);
                              
               char j_line_str[LINE_BYTES];
               int  j_bytes_to_read = LINE_BYTES; // Arbitrary... Needs to be longer than the longest line in the file.
               while(fgets (j_line_str, j_bytes_to_read, fj) != NULL)             
               {
                  // Inner line specific variables here
                  char j_position[POSITION_STR_LEN];
                  int j_umi_count     = 0;
                  int j_token_counter = 1;
                  char * jtoken       = NULL;   
              
                  // Eliminate newline
                  char *newline = strchr( j_line_str, '\n' );
                  if ( newline )
                  *newline = 0;
 
                  // This gets tokens of the inner string
                  char * jptr;
                  jtoken = strtok_r(j_line_str, sep, &jptr);                     

                  while( jtoken != NULL )
                  {                         
                     if(j_token_counter == 1)
                     {
                        strcpy(j_position,jtoken);
                     }
                     else if(j_token_counter == 2)
                     {
                        j_umi_count = atoi(jtoken);
                     }
                     else
                     {
                        // j_token_counter must be > 2 and therefore jtoken is a UMI
                        //printf( "\t\tJ: %d   %s   %d   %s\n", j_token_counter,j_position,j_umi_count,jtoken );
 
                        // Is the i UMI the same as the j UMI? If we get a match, then we add the j position to the
                        // list of positions belonging to the i UMI

                        if(!(strcmp (itoken,jtoken)))
                        {
                           //printf( "\t\tMATCH: %s   %s   %s   %s\n", i_position,itoken,j_position,jtoken );
                           
                           //
                           // MAIN OUTPUT
                           //
                           printf("%s ",j_position);
                        }
                     }
 
                     jtoken = strtok_r(NULL, sep, &jptr);
                     j_token_counter++;
 
                  } // Inner loop while jtoken
                                     
                  j_line_counter++;
                  
                } // Inner loop while getline
                
                //
                // MAIN OUTPUT
                //
                printf("\n");
                                                                           
            } // else it is a new i UMI
   
         } // else itoken is a UMI
                  
         i_token_counter++;
         itoken = strtok_r(NULL, sep, &iptr);
         
      } // Outer loop while itoken
      
      // This is the end of the outer line. Reset.
      i_token_counter = 1;
      i_line_counter++;
      
   } // Outer loop while getline
   
   // Sort the keys
   HASH_SORT(umi_counts, string_sort);
   
   // Logging output
   FILE *flog = fopen("analyzer_umis.txt", "w");
   
   // Iterate the sorted hash and print the UMI
   struct count_strings *s;

   for(s=umi_counts; s != NULL; s=s->hh.next) 
   {
      fprintf(flog,"%s\n", s->umi);
   }   
      
   fprintf(stderr,"Closing files\n");
   
   fclose(flog); 
   fclose(fi);
   fclose(fj);
  
   return 0;
}


// Comparison function for strings, used by HASH_SORT
int string_sort(struct count_strings *a, struct count_strings *b) 
{
    return strcmp(a->umi,b->umi);
}


// Get the number of lines in a file
int get_num_lines(char * file)
{
   char str[100];
   sprintf(str,"wc -l %s", file);
   FILE *fw = popen (str,"r");

   if (fw == NULL)
   {
        fprintf(stderr, "Unable to open %s for line counting\n", file);
        return(-1);
   }

   char wc_output[100];
   fgets (wc_output,sizeof(wc_output),fw);
   char * num_lines_str = NULL;
   num_lines_str = strtok(wc_output, " ");
   int num_lines = atoi(num_lines_str);
   fclose(fw);
   return num_lines;
}











