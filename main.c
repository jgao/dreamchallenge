#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int KI_LENGTH = 15; // first 15 of 31 in 31-mer
int KR_LENGTH = 16; // last 16 of 31 in 31-mer

int KMER_LENGTH = 31;


// from Pollux
void deleteCharacter(char* string, int pos)
{
    int length = strlen(string);
    
    if (pos < 0 || pos >= length)
        return;
    
    for(int i = pos; i < length; i++)
    {
        string[i] = string[i + 1];
    }
    
    string[length] = '\0';
}
void trimSpaces(char* string)
{
    // Delete leading spaces.
    while(strlen(string) >= 1 && isspace(string[0]))
    {
        deleteCharacter(string, 0);
    }
    
    // Delete trailing N's.
    while(strlen(string) >= 1 && isspace(string[strlen(string) - 1]))
    {
        deleteCharacter(string, strlen(string) - 1);
    }
}

void printBinary(unsigned long long int n){

	int i = 0;

    while (i<64) {
	    if (n & 1)
	        printf("1");
	    else
	        printf("0");

	    if (i%4 == 3) printf(" ");

	    n >>= 1;
	    i++;
	}
	printf("\n");

}

unsigned int nucleotideAlphaToNumeric(char nucleotide) {
	if (nucleotide == 'A') return 0x0;
	if (nucleotide == 'G') return 0x1;
	if (nucleotide == 'C') return 0x2;
	if (nucleotide == 'T') return 0x3;

	printf("error: could not match nucleotide %c\n", nucleotide);
	return 0x0;
}

unsigned int nucleotideLexiRank(unsigned long long int nucleotide) {
	if (nucleotide == 0x0000000000000000) return 1;
	if (nucleotide == 0x4000000000000000) return 3;
	if (nucleotide == 0x8000000000000000) return 2;
	if (nucleotide == 0xC000000000000000) return 4;

	printf("error: could not match nucleotide %d\n", nucleotide);
	return 0;
}

unsigned long long int getKMerNumeric(char* kmer_alpha)
{
    unsigned long long int result = 0x0;
    
	for (int i = 0; i < KMER_LENGTH; i++) {
		result = (result | nucleotideAlphaToNumeric(kmer_alpha[i])) << 2;
	}

	return result;
}


/*
unsigned int getKI(unsigned long long int kmer)
{
    const unsigned long long int KIMASK = 0x00000000FFFFFFFF;
    
    return (kmer & KIMASK) >> 2; // we get rid of the last 2 bits 
}

unsigned int getKR(unsigned long long int kmer)
{
    const unsigned long long int KRMASK = 0xFFFFFFFF00000000;

    return (kmer & KRMASK) >> 32; // need to shift over so casting doesn't truncate the part we want
}
*/


const unsigned long long int HIMASK = 0xFFFFFFFF00000000;
const unsigned long long int LOMASK = 0x00000000FFFFFFFF;

// Note M2 is the one with the 00 ending
const unsigned long long int M1 = 0x961C69E3;
const unsigned long long int M2 = 0x69E3961C;

unsigned int getKI(unsigned long long int kmer)
{   
	unsigned int hi = (kmer & HIMASK) >> 32;
	unsigned int lo = (kmer & LOMASK);

    return ((hi & M2) + (lo & M1)) >> 2; // we get rid of the last 2 bits
}

unsigned int getKR(unsigned long long int kmer)
{
	unsigned int hi = (kmer & HIMASK) >> 32;
	unsigned int lo = (kmer & LOMASK);

    return ((hi & M1) + (lo & M2)); 
}

// returns -1 if k1 is smaller OR EQUAL, 1 if k2 is smaller
unsigned int getLexiSmallerKMer(unsigned long long int k1, unsigned long long int k2)
{
	const unsigned long long int MASK = 0xC000000000000000;

	for (int i = 0; i < KMER_LENGTH; i++) {
		if (nucleotideLexiRank(k1 & MASK) < nucleotideLexiRank(k2 & MASK)) return -1;
		if (nucleotideLexiRank(k1 & MASK) > nucleotideLexiRank(k2 & MASK)) return 1;

		k1 = k1 << 2;
		k2 = k2 << 2;
	}
	// everything is the same
	return 0;
}

void writeAsNucleotides(FILE* file, unsigned long long int sequence)
{
    const unsigned long long int MASK = 0x3;

    // Iterate over all nucleotides:
    for(int i = 0; i < 31; i++)
    {
        unsigned int position = (i * 2) % 64;     // bit position w/i sequence
        
        unsigned long long int nucleotide = (sequence >> (62 - position)) & MASK;
        
        // Print appropriate character:
        if(nucleotide == 0x0)
        {
            fprintf(file, "A");
        }
        else if(nucleotide == 0x1)
        {
            fprintf(file, "G");
        }
        else if(nucleotide == 0x2)
        {
            fprintf(file, "C");
        }
        else if(nucleotide == 0x3)
        {
            fprintf(file, "T");
        }
        // Error:
        else
        {
            fprintf(file, "E");
        }
    }
}
void printAsNucleotides(unsigned long long int sequence)
{
    writeAsNucleotides(stdout, sequence);
}

unsigned long long int getReverse(unsigned long long int sequence)
{
    unsigned long long int MASK = 0x3;
    
    unsigned long long int nucleotide;
    unsigned long long int reverse = 0x0;

    
    for(int i = 0; i < 32; i++)
    {
        // Get the last nucleotide:
        nucleotide = (sequence >> (i * 2)) & MASK;        
                
        // Add the nucleotide to the reverse:
        reverse = reverse | (nucleotide << ((32 - i - 1) * 2));
    }  
    
    return reverse;
}

unsigned long long int createReverseCompliment(unsigned long long int kmer)
{
    
    unsigned long long int reverse;
    
    unsigned long long int reverseCompliment;

    reverse = getReverse(kmer);
    reverseCompliment = ~(reverse);
    return reverseCompliment << 2;
}


int main() {
	FILE* dump_file = fopen("dump", "r");

	unsigned int num_possible_ki = 1 << (KI_LENGTH * 2); // shortcut to 2 ^ (KI_LENGTH * 2)

	printf("%d\n", num_possible_ki);

	unsigned long long int* ki_index = (unsigned long long int*)malloc(num_possible_ki * sizeof(unsigned long long int));
	// this counter will be used to remember the next spot to insert a KR
	unsigned long long int* temp_kr_list_bookmark = (unsigned long long int*)malloc(num_possible_ki * sizeof(unsigned long long int));
	for (unsigned int i = 0; i < num_possible_ki; i++) {
		ki_index[i] = 0;
		temp_kr_list_bookmark[i] = 0;
	}

	// first pass of file: read occurences of each KI
	char count_line[KMER_LENGTH + 10];
	char kmer_alpha[KMER_LENGTH + 10];
	int ccc = 0;
	while (fgets(count_line, KMER_LENGTH + 10, dump_file) != NULL
		&& fgets(kmer_alpha, KMER_LENGTH + 10, dump_file) != NULL) {
		
		trimSpaces(count_line);
		deleteCharacter(count_line, 0);

		trimSpaces(kmer_alpha);

		unsigned long long int kmer_numeric = getKMerNumeric(kmer_alpha);

		unsigned int ki = getKI(kmer_numeric);
		//unsigned int kr = getKR(kmer_numeric);

		ki_index[ki]++;

		// short iterations for testing
		//if (ccc++ == 5) break;
	}

	fclose(dump_file);

	unsigned long long int cumulative_sum = 0;
	// we want a cumulative index, so we do a roll-up sum
	for (unsigned int i = 0; i < num_possible_ki; i++) {
		// the i-1th index should already be rolled up
		unsigned long long int temp = ki_index[i];
		ki_index[i] = cumulative_sum;
		cumulative_sum = cumulative_sum + temp;
	}
printf("===");
	// the final roll up is just the total number of k-mers present
	unsigned long long int num_kmer = cumulative_sum;

	// second pass we build the counters by each KR
	dump_file = fopen("dump", "r");
	
	unsigned int* kr_list = (unsigned int*)malloc(num_kmer * sizeof(unsigned int));
	unsigned int* kmer_frequency = (unsigned int*)malloc(num_kmer * sizeof(unsigned int));
	// TODO: add another counter for dream challenge stuff
ccc = 0;
	while (fgets(count_line, KMER_LENGTH + 10, dump_file) != NULL
		&& fgets(kmer_alpha, KMER_LENGTH + 10, dump_file) != NULL) {
		
		trimSpaces(count_line);
		deleteCharacter(count_line, 0);

		trimSpaces(kmer_alpha);

		int kmer_count = (int)strtod(count_line, NULL);

		unsigned long long int kmer_numeric = getKMerNumeric(kmer_alpha);

		unsigned int ki = getKI(kmer_numeric);
		unsigned int kr = getKR(kmer_numeric);

		unsigned long long int kr_list_start_pos = ki_index[ki];
		
		kr_list[kr_list_start_pos + temp_kr_list_bookmark[ki]] = kr;
		kmer_frequency[kr_list_start_pos + temp_kr_list_bookmark[ki]] = kmer_count;

		temp_kr_list_bookmark[ki]++;

		// short iterations for testing
		//if (ccc++ == 5) break;
	}

    /*FILE* histo_file = fopen("jellyhash_bucket_histo", "w");
	for (unsigned int i = 0; i < num_possible_ki; i++) {
        fprintf(histo_file, "%llu %llu\n", i, temp_kr_list_bookmark[i]);
        if (i % 25000000 == 0) printf("%d / %d \n", i, num_possible_ki);
    }*/

printf("ready\n");
	while(1) {
		char search[32];

		scanf("%s", search);

		if(search[0] == '!') break;

		unsigned long long int kmer_numeric = getKMerNumeric(search);
		unsigned long long int rev_comp = createReverseCompliment(kmer_numeric);
		printAsNucleotides(kmer_numeric);printf("kmer\n");
		printAsNucleotides(rev_comp);printf("revcomp\n");

		if (getLexiSmallerKMer(kmer_numeric, rev_comp) == 1)
			kmer_numeric = rev_comp;

		printf("This is the lexi_smaller one: ");printAsNucleotides(kmer_numeric);printf("\n");
		unsigned int ki = getKI(kmer_numeric);
		unsigned int kr = getKR(kmer_numeric);
		printf("%llu ", ki); printBinary(ki);
		printf("%llu ", kr);printBinary(kr);

		unsigned long long int kr_list_start_pos = ki_index[ki];
		int is_found = 0;

		printf("%llu\n", kr_list_start_pos);
		for (unsigned int i = 0; i < temp_kr_list_bookmark[ki]; i++) {
			if (kr_list[kr_list_start_pos + i] == kr) {
				printf("%llu\n", kmer_frequency[kr_list_start_pos + i]);
				is_found = 1;
				break;
			}
		}

		if (!is_found) {
			printf("KMER NOT FOUND\n");
		}

	}

	fclose(dump_file);

	free(ki_index);
	free(temp_kr_list_bookmark);
	free(kr_list);
	free(kmer_frequency);

	return 0;
}