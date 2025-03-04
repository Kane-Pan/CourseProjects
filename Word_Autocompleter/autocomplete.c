#include "autocomplete.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int compare_terms(const void *a, const void *b) 
{
    const term *termA = (const term *)a;
    const term *termB = (const term *)b;
    return strcmp(termA->term, termB->term);
}

int compare_weights(const void *a, const void *b)
{
    const term *termA = (const term *)a;
    const term *termB = (const term *)b;
    if (termA->weight > termB->weight)
        return -1;
    else if (termA->weight < termB->weight)
        return 1;
    else
        return 0;
}

void read_in_terms(struct term **terms, int *pnterms, char *filename)
{
    FILE *fp = fopen(filename, "r");
    int n;
    char line[256]; // assume the maximum characters per line is 256

    fscanf(fp, "%d", &n); //get the number of characters n from the first number in txt
    fgetc(fp); // moves reader to the first line starting where cities are

    *terms = malloc(n * sizeof(term));
    for (int i = 0; i < n; i++) {
        fgets(line, sizeof(line), fp);
        // scan the line then place the respective contents into weight and term in terms
        sscanf(line, "%lf %[^\n]", &((*terms)[i].weight), (*terms)[i].term);
    }
    fclose(fp);

    // sort lexicographically
    qsort(*terms, n, sizeof(term), compare_terms);
    *pnterms = n;
}

int lowest_match(term *terms, int nterms, char *substr)
{
    // binary search
    int left = 0;
    int right = nterms - 1;

    while (left <= right) {
        int mid = (left + right) / 2;
        int cmp = strncmp(terms[mid].term, substr, strlen(substr));

        if (cmp == 0) {
            while (mid > 0 && strncmp(terms[mid - 1].term, substr, strlen(substr)) == 0){
                // move down leftward to find least
                mid--;
            }
            return mid;
        }
        else if (cmp > 0){
            right = mid - 1;
        }
        else {
            left = mid + 1;
        }
    }
    return -1;
}

int highest_match(term *terms, int nterms, char *substr)
{
    int left = 0;
    int right = nterms - 1;

    while (left <= right) {
        int mid = (left + right) / 2;
        int cmp = strncmp(terms[mid].term, substr, strlen(substr));

        if (cmp == 0) {
            while (mid < nterms - 1 && strncmp(terms[mid + 1].term, substr, strlen(substr)) == 0){
                // move right to find lexicographical highest
                mid++;
            }
            return mid;
        }
        else if (cmp > 0){
            right = mid - 1;
        }
        else {
            left = mid + 1;
        }
    }
    return -1;
}

void autocomplete(term **answer, int *n_answer, term *terms, int nterms, char *substr)
{
    int highest_index = highest_match(terms, nterms, substr);
    int lowest_index = lowest_match(terms, nterms, substr);

    // since it is lexicographically sorted, the words in between will also have the substr
    if (lowest_index == -1 || highest_index == -1) {
        *n_answer = 0;
        return;
    }

    *n_answer = highest_index - lowest_index + 1;
    *answer = malloc((*n_answer) * sizeof(term));

    for (int i = 0; i < *n_answer; i++){
        (*answer)[i] = terms[lowest_index + i];
    }

    // sort by weights
    qsort(*answer, *n_answer, sizeof(term), compare_weights);
}
