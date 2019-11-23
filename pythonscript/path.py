def add(x, y):
    try:
        from nltk.corpus import wordnet
        return wordnet.synsets(x)[0].path_similarity(wordnet.synsets(y)[0])
    except IndexError:
        return 0.0