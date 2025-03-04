#include "autocomplete.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
    term *terms = NULL;
    int nterms = 0;
    char *filename = "cities.txt";

    // Read in all the terms from the file.
    read_in_terms(&terms, &nterms, filename);
    printf("Read %d terms from %s\n", nterms, filename);

    // Prompt for a query string.
    char query[200];
    printf("Enter query string: ");
    if (fgets(query, sizeof(query), stdin) != NULL) {
        query[strcspn(query, "\n")] = '\0';  // Remove trailing newline.
    }

    // Call autocomplete with the given query.
    term *answers = NULL;
    int n_answers = 0;
    autocomplete(&answers, &n_answers, terms, nterms, query);

    // Print the results.
    if (n_answers == 0) {
        printf("No matches found for \"%s\".\n", query);
    } else {
        printf("Found %d matches:\n", n_answers);
        for (int i = 0; i < n_answers; i++) {
            printf("%s (weight: %.2lf)\n", answers[i].term, answers[i].weight);
        }
    }

    // Free allocated memory.
    free(terms);
    free(answers);

    getchar();
    return 0;
}
