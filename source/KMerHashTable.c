/*

Pollux
Copyright (C) 2014  Eric Marinier

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <stdlib.h>
#include "KMerHashTable.h"
#include "Utility.h"
#include "Correction.h"

const int KI_LENGTH = 15; // first 15 of 31 in 31-mer
const int KR_LENGTH = 16; // last 16 of 31 in 31-mer

const int KMER_LENGTH = 31;

const unsigned long long int HIMASK = 0xFFFFFFFF00000000;
const unsigned long long int LOMASK = 0x00000000FFFFFFFF;

// Note M2 is the one with the 00 ending
const unsigned long long int M1 = 0x55555553;
const unsigned long long int M2 = 0xAAAAAAAC;

static inline void printProgress(int x, int n, int r)
{
    // Only update r times.
    if (n/r == 0 || x % (n/r) != 0) return; // Mod 0 causes problems.
    
    double result;
    
    if(n / r != 1)
        result = x / (n/r) * (100/r);
    else
        result = (x + 1) / ((double)n/(double)r) * ((double)100/(double)r);
 
    printf("%d%% ", (int)result);
    fflush(stdout);
}

unsigned long long int getKMerNumeric(char* kmer_alpha)
{
    unsigned long long int result = 0x0;
    
    for (int i = 0; i < KMER_LENGTH; i++) {
        result = (result | nucleotideAlphaToNumeric(kmer_alpha[i])) << 2;
    }

    return result;
}

unsigned int getKI(unsigned long long int kmer)
{
    // const unsigned long long int KIMASK = 0x00000000FFFFFFFF;
    
    unsigned int hi = (kmer & HIMASK) >> 32;
    unsigned int lo = (kmer & LOMASK);

    return ((hi & M2) + (lo & M1)) >> 2; // we get rid of the last 2 bits
}

unsigned int getKR(unsigned long long int kmer)
{
    // const unsigned long long int KRMASK = 0xFFFFFFFF00000000;

    unsigned int hi = (kmer & HIMASK) >> 32;
    unsigned int lo = (kmer & LOMASK);

    return ((hi & M1) + (lo & M2)); // we get rid of the last 2 bits 
}

unsigned int getlowKMerThreshold(KMerHashTable* kmerTable) {

    unsigned long long int MAX_KMER_COUNT = (1024 + 1);
    unsigned int* dist = createDistribution(kmerTable, MAX_KMER_COUNT);

    // We know the kmer count for 0 and 1 are both 0, so start at 2
    unsigned int currentKMerCount = 2;
        
    // Loop until we find a low-to-high number transitions:
    while (currentKMerCount <= MAX_KMER_COUNT && dist[currentKMerCount] > dist[currentKMerCount + 1])
    {
        currentKMerCount++;
    }

    free(dist);

    return currentKMerCount;
}

KMerHashTable* newKMerHashTable(char* jellyfishFile)
{
    // only need to save half
    uint32_t sizeof_revcomp_dict = 0xFFFFFFFF;
    uint32_t* revcomp_dict = (uint32_t*)malloc(sizeof(uint32_t) * sizeof_revcomp_dict);
    
    for (uint32_t i = 0; i < sizeof_revcomp_dict; i++) {
        revcomp_dict[i] = ~getReverse16(i);
    }




printf("START OF BUILDING JELLYFISH HASHTABLE\n");
    KMerHashTable* table = (KMerHashTable*)malloc(sizeof(KMerHashTable));
 
    table->num_possible_ki = 1 << (KI_LENGTH * 2); // shortcut to 2 ^ (KI_LENGTH * 2)
    
    table->ki_index = (unsigned long long int*)malloc(table->num_possible_ki * sizeof(unsigned long long int));
    table->num_kr_per_ki = (unsigned long long int*)malloc(table->num_possible_ki * sizeof(unsigned long long int));

    // TODO: add another counter for dream challenge stuff

    FILE* jellyfish_file = fopen(jellyfishFile, "r");

    for (int i = 0; i < table->num_possible_ki; i++) {
        table->ki_index[i] = 0;
        table->num_kr_per_ki[i] = 0;
    }

    // first pass of file: read occurences of each KI
    char count_line[KMER_LENGTH + 10];
    char kmer_alpha[KMER_LENGTH + 10];

int iter=0;
    while (fgets(count_line, KMER_LENGTH + 10, jellyfish_file) != NULL
        && fgets(kmer_alpha, KMER_LENGTH + 10, jellyfish_file) != NULL) {
        
        trimSpaces(count_line);
        deleteCharacter(count_line, 0);

        trimSpaces(kmer_alpha);

        int kmer_count = (int)strtod(count_line, NULL);

        // ignore all kmers appearing once
    //if (kmer_count == 1) continue;

        unsigned long long int kmer_numeric = getKMerNumeric(kmer_alpha);
//unsigned long long int reverse_comp = createReverseComplimentKMer(kmer_numeric);

        unsigned int ki = getKI(kmer_numeric);
//unsigned int rcki = getKI(reverse_comp);

        table->ki_index[ki]++;
//table->ki_index[rcki]++;

        /*printBinary(kmer_numeric);
        printf("%d / %d\n", (0xFFFFFFFF00000000 & kmer_numeric) >> 32, sizeof_revcomp_dict);
        printBinary(revcomp_dict[(0xFFFFFFFF00000000 & kmer_numeric) >> 32]);
        printf("%d / %d\n", 0x00000000FFFFFFFF & kmer_numeric, sizeof_revcomp_dict);
        printBinary(revcomp_dict[0x00000000FFFFFFFF & kmer_numeric]);
*/
        unsigned long long int rc_hi = revcomp_dict[(0xFFFFFFFF00000000 & kmer_numeric) >> 32];
        unsigned long long int rc_lo = revcomp_dict[0x00000000FFFFFFFF & kmer_numeric];

        printBinary((rc_hi << 2));
        printBinary((rc_lo << 34));

        unsigned long long int rc = (rc_hi << 2) | (rc_lo << 34);

        printBinary(rc);
        printBinary(createReverseComplimentKMer(kmer_numeric));
        //printBinary(getKMerNumeric("TAAATTGACCATCAACCGCACCCGACATCAC"));
        printf("\n");



        if (iter++ >= 4) break;
    }
exit(0);

    fclose(jellyfish_file);
printf("DONE FIRST PASS\n");
    unsigned long long int cumulative_sum = 0;
    // we want a cumulative index, so we do a roll-up sum
    for (unsigned int i = 0; i < table->num_possible_ki; i++) {
        // the i-1th index should already be rolled up
        unsigned long long int temp = table->ki_index[i];
        table->ki_index[i] = cumulative_sum;
        cumulative_sum = cumulative_sum + temp;
    }

    // total amount of kmers is just the final cumulative_sum
    table->num_kmer = cumulative_sum;

    table->kr_list = (unsigned int*)malloc(table->num_kmer * sizeof(unsigned int));
    table->kmer_frequency = (unsigned int*)malloc(table->num_kmer * sizeof(unsigned int));

    // second pass we build the counters by each KR
    jellyfish_file = fopen(jellyfishFile, "r");

    while (fgets(count_line, KMER_LENGTH + 10, jellyfish_file) != NULL
        && fgets(kmer_alpha, KMER_LENGTH + 10, jellyfish_file) != NULL) {
        
        trimSpaces(count_line);
        deleteCharacter(count_line, 0);

        trimSpaces(kmer_alpha);

        int kmer_count = (int)strtod(count_line, NULL);

        // ignore all kmers appearing once
        if (kmer_count == 1) continue;

        unsigned long long int kmer_numeric = getKMerNumeric(kmer_alpha);
//unsigned long long int reverse_comp = createReverseComplimentKMer(kmer_numeric);

        unsigned int ki = getKI(kmer_numeric);
        unsigned int kr = getKR(kmer_numeric);
//unsigned int rcki = getKI(reverse_comp);
//unsigned int rckr = getKR(reverse_comp);

        unsigned long long int kr_list_start_pos = table->ki_index[ki];
//unsigned long long int rckr_list_start_pos = table->ki_index[rcki];
        
        table->kr_list[kr_list_start_pos + table->num_kr_per_ki[ki]] = kr;
        table->kmer_frequency[kr_list_start_pos + table->num_kr_per_ki[ki]] = kmer_count;
//table->kr_list[rckr_list_start_pos + table->num_kr_per_ki[rcki]] = rckr;
//table->kmer_frequency[rckr_list_start_pos + table->num_kr_per_ki[rcki]] = kmer_count;

        table->num_kr_per_ki[ki]++;
//table->num_kr_per_ki[rcki]++;
    }

    fclose(jellyfish_file);
printf("END OF BUILDING JELLYFISH HASHTABLE\n");

    return table;
}

unsigned int nucleotideLexiRank(unsigned long long int nucleotide) {
    if (nucleotide == 0x0000000000000000) return 1;
    if (nucleotide == 0x4000000000000000) return 3;
    if (nucleotide == 0x8000000000000000) return 2;
    if (nucleotide == 0xC000000000000000) return 4;

    printf("error: could not match nucleotide %d\n", nucleotide);
    return 0;
}

// returns -1 if k1 is smaller OR EQUAL, 1 if k2 is smaller
unsigned int getLexiSmallerKMer(unsigned long long int k1, unsigned long long int k2)
{
    const unsigned long long int MASK = 0xC000000000000000;

    for (int i = 0; i < 31; i++) {
        if (nucleotideLexiRank(k1 & MASK) < nucleotideLexiRank(k2 & MASK)) return -1;
        if (nucleotideLexiRank(k1 & MASK) > nucleotideLexiRank(k2 & MASK)) return 1;

        k1 = k1 << 2;
        k2 = k2 << 2;
    }
    // everything is the same
    return 0;
}

// Helper function that does the actual lookup, without reverse compliment checks.
unsigned long long int KMerTableLookupSingleDirection(KMerHashTable* kmerTable, unsigned long long int kmer)
{
    unsigned int ki = getKI(kmer);
    unsigned int kr = getKR(kmer);

    unsigned long long int kr_list_start_pos = kmerTable->ki_index[ki];
    int count = 0;

    for (unsigned int i = 0; i < kmerTable->num_kr_per_ki[ki]; i++) {
        if (kmerTable->kr_list[kr_list_start_pos + i] == kr) {
            count = kmerTable->kmer_frequency[kr_list_start_pos + i];
            break;
        }
    }
    
    // because we discard unique kmers, we (optimistically) assume that uncatalogued kmers have a count of 1
    if(count == 0)
    {
        count = 1;
    }
    
    return count;
}

unsigned long long int KMerTableLookup(KMerHashTable* kmerTable, unsigned long long int kmer)
{
    /**/unsigned long long int reverse_comp = createReverseComplimentKMer(kmer);

    if (getLexiSmallerKMer(kmer, reverse_comp) == 1)
        kmer = reverse_comp;
/**/
    unsigned long long int count = KMerTableLookupSingleDirection(kmerTable, kmer);
    return count;
}

unsigned int getMaxKMerCount(KMerHashTable* kmerTable)
{
    unsigned int max = -1;
    printf("inside get maxcount %llu\n", kmerTable->num_kmer);
    for (unsigned long long int i = 0; i < kmerTable->num_kmer; i++) {
        max = getMax(max, kmerTable->kmer_frequency[i]);
    }
    
    return max;
}

unsigned int* createDistribution(KMerHashTable* kmerTable, unsigned int max)
{    
    unsigned int* distribution = (unsigned int*)malloc((max + 1) * sizeof(unsigned int));
    
    // Initialize:
    for(int i = 0; i <= max; i++)
    {
        distribution[i] = 0;
    }

    for (unsigned long long int i = 0; i < kmerTable->num_kmer; i++) {
        if (kmerTable->kmer_frequency[i] > max) continue;

        distribution[kmerTable->kmer_frequency[i]]++;
    }
    
    return distribution;
}

unsigned int getNumRepeats(KMerHashTable* kmerTable)
{   
    // Variables:
    unsigned int repeats = 0;
    unsigned int max = getMaxKMerCount(kmerTable);
    
    // Data Structures:
    unsigned int* distribution = createDistribution(kmerTable, max);    

    for(int i = 2; i <= max; i++)
    {
        repeats += distribution[i] * i;
    }
    
    free(distribution);
    
    return repeats;
}


void freeKMerHashTable(KMerHashTable* kmerTable) {
    free(kmerTable->ki_index);
    free(kmerTable->num_kr_per_ki);
    free(kmerTable->kr_list);
    free(kmerTable->kmer_frequency);

    free(kmerTable);
}
