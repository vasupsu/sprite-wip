#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<ctype.h>
#include "faiTrie.h"

//TrieNode *trie = NULL;

TrieNode *newTrieNode (char ch, int len, int contigNo)
{
   TrieNode *node = (TrieNode *)malloc(sizeof(TrieNode));
   node->ch=ch;
   node->contigNo=contigNo;
   node->length=len;
   node->numMatch=0;
   node->next = (TrieNode **)malloc(36 * sizeof(TrieNode *));
   bzero (node->next, 36 * sizeof(TrieNode *));
   return node;
}

void insertStr (TrieNode *trie, char *str, int len, int cNo)
{
   int i;
   
   TrieNode *iter = trie;
   for (i=0; str[i]!='\0'; i++)
   {
      if (!isalnum(str[i])) continue;
      char ch=str[i];
      int index=-1;
      if (isalpha(ch))
         if (islower(ch))
            index=ch-'a'+10;
         else
            index=ch-'A'+10;
      else
         index=ch-'0';
      iter->numMatch++;
      if (iter->next[index] == NULL)
      {
         iter->next[index]=newTrieNode(ch, -1, cNo);
      }
      iter = iter->next[index];
   }
   iter->length=len;
   iter->numMatch++;
}

char trieStr[100];
int strIndex=0;
void printTrie (TrieNode *head)
{
   trieStr[strIndex++]=head->ch;
   if (head->length!=-1)
   {
      trieStr[strIndex]='\0';
      printf ("%s\t%d\n", trieStr, head->length);
   }
   int i=0;
   // int found=0;
   for(i=0; i<36; i++)
   {
      if (head->next[i]!=NULL)
      {
         printTrie(head->next[i]);
         // found=1;
      }
   }
   strIndex--;
}

int getContigNo(TrieNode *head, char *str, int index)
{
   if (str[index]=='\0')
      return head->contigNo;
   if (!isalnum(str[index])) 
      return getContigNo(head, str, index+1);
   char ch=str[index];
   int nIndex=-1;
   if (isalpha(ch))
      if (islower(ch))
         nIndex=ch-'a'+10;
      else
         nIndex=ch-'A'+10;
   else
      nIndex=ch-'0';
   return getContigNo(head->next[nIndex], str, index+1);
}

/*int main()
{
   trie = newTrieNode(' ', -1, NULL);
   char contigName[100];
   int length;
   char line[100];
   char faiFile[]="/gpfs/home/vxr162/scratch/hashdata/ucsc.hg19.fasta.fai";
   FILE *fp = fopen(faiFile, "r");
   fgets(line, 100, fp);
   while (!feof(fp))
   {
      sscanf(line, "%s\t%d", contigName, &length);
//      printf ("*%s*%d*\n", contigName, length);
      insertStr (contigName, length, NULL);
      fgets(line, 100, fp);
   }
   fclose(fp);
   printTrie (trie);
   int len=getLength(trie, "chrUn_gl000236", 0);
   printf ("chrUn_gl000236::%d\n", len);
}*/
