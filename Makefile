
CXXFLAGS = -std=c++17 -Wextra -Wall -Werror -pedantic -Wunused -Wshadow

.PHONY: run clean zip

pro: pro.cpp
	mpic++ -o pro pro.cpp $(CXXFLAGS)

clean:
	rm -f pro xkrejc69.zip

zip: pro
	zip xkrejc69.zip pro.cpp test.sh xkrejc69.pdf
