import math
import os
os.chdir(r"C:\Users\kaney\OneDrive\Documents\University\ESC180\Final Project Synonyms")

def norm(vec):
    '''Return the norm of a vector stored as a dictionary, as
    described in the handout for Project 3.
    '''

    sum_of_squares = 0.0
    for x in vec:
        sum_of_squares += vec[x] * vec[x]

    return math.sqrt(sum_of_squares)


def cosine_similarity(vec1, vec2):
    similarity = 0.0
    nums1 = list(vec1.values())
    nums2 = list(vec2.values())
    nums1squared = []
    nums2squared = []

    # Find the magnitude of the two vectors
    for i in range(len(nums1)):
        nums1squared.append(nums1[i] ** 2)
    for i in range(len(nums2)):
        nums2squared.append(nums2[i] ** 2)

    vec1mag = math.sqrt(sum(nums1squared))
    vec2mag = math.sqrt(sum(nums2squared))

    # Avoid division by zero
    if vec1mag == 0 or vec2mag == 0:
        return -1

    # Find dot product
    dotproduct = 0
    for word, freq in vec1.items():
        if word in vec2.keys():
            dotproduct += vec1[word] * vec2[word]

    similarity = dotproduct / (vec1mag * vec2mag)

    return similarity


def build_semantic_descriptors(sentences):
    # Dictionary to hold the semantic descriptors
    semantic_descriptors = {}

    for sentence in sentences:
    # Set to avoid duplicate words
        unique_words = set(sentence)

        # Update co-occurrence counts for each word in the sentence
        for word in unique_words:
            word = word.lower()
            if word not in semantic_descriptors:
                semantic_descriptors[word] = {}
            for sep_word in unique_words:
                if word != sep_word:
                    if sep_word not in semantic_descriptors[word]:
                        semantic_descriptors[word][sep_word] = 0
                    semantic_descriptors[word][sep_word] += 1

    return semantic_descriptors


def build_semantic_descriptors_from_files(filenames):
    sentences = []

    for i in range(len(filenames)):
        current_file = open(filenames[i], "r", encoding="latin1")
        content = current_file.read().lower()

        # Replace sentence-ending punctuation with a unique marker
        for sep in [".", "!", "?"]:
            content = content.replace(sep, "|")

        # Split the text into sentences
        raw_sentences = content.split("|")

        for raw_sentence in raw_sentences:
            # Split the sentence into words
            words = raw_sentence.split()

            # Remove unwanted punctuation from the words
            cleaned_words = []
            for word in words:
                cleaned_words.append(word.strip(",-:;"))

            if cleaned_words:
                sentences.append(cleaned_words)

        current_file.close()

    return build_semantic_descriptors(sentences)



def most_similar_word(word, choices, semantic_descriptors, similarity_fn):
    max_similarity_score = -1
    list(filter(lambda a: a != "", choices))
    best_choice = choices[0]

    if word not in semantic_descriptors.keys():
        return best_choice

    sem_desc_of_word = semantic_descriptors[word]

    for choice in choices:
        # If word is not in the dictionary
        if choice not in semantic_descriptors.keys():
            continue

        sem_desc_of_choice = semantic_descriptors[choice]
        sim_score_of_choice = similarity_fn(sem_desc_of_word, sem_desc_of_choice)

        if sim_score_of_choice > max_similarity_score:
            max_similarity_score = sim_score_of_choice
            best_choice = choice

    return best_choice


def run_similarity_test(filename, semantic_descriptors, similarity_fn):

    test_file = open(filename, "r", encoding="latin1")
    line_number = 1
    correct = 0
    total_num_of_questions = 0

    for question in test_file:
        word_position = 1
        words = question.split()  # Split the question into words
        choices = []

        for word in words:
            if word_position == 1:
                question_word = word
            elif word_position == 2:
                answer = word
            else:
                choices.append(word)
            word_position += 1  # Increment word position

        if most_similar_word(question_word, choices, semantic_descriptors, similarity_fn) == answer:
            correct += 1

        total_num_of_questions += 1

    test_file.close()

    return (correct / total_num_of_questions) * 100


if __name__ == "__main__":



    sem_descriptors = build_semantic_descriptors_from_files(["reference_swannsway.txt", "reference_warandpeace.txt", "reference_prideandprejudice.txt"])
    res = run_similarity_test("testlist.txt", sem_descriptors, cosine_similarity)
    print(res, "percent of the guesses were correct")










