typedef struct TrieNode
{
   char ch;
   int contigNo;//contig posTable
   int length;//contig length
   int numMatch;//# of strings with common prefix from root to this node
   struct TrieNode **next;
}TrieNode;

extern TrieNode *newTrieNode (char ch, int len, int contigNo);
extern void insertStr (TrieNode *trie, char *str, int len, int cNo);
extern int getContigNo(TrieNode *head, char *str, int index);
