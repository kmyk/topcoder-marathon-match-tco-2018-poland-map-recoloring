.PHONY: build
.DEFAULT: build

PROBLEM := MapRecoloring
CXX := g++
CXXFLAGS := -std=c++11 -Wall -O2 -g -DLOCAL

build: a.out tester.jar

run: a.out tester.jar
	java -jar tester.jar -exec ./a.out

a.out: main.cpp ${PROBLEM}.cpp
	${CXX} ${CXXFLAGS} $<

tester.jar: ${PROBLEM}Vis.java
	javac $<
	jar cvfe tester.jar ${PROBLEM}Vis *.class

URL := https://community.topcoder.com/longcontest/?module=ViewProblemStatement&rd=17149&pm=14893
submit:
	oj submit '${URL}' --language C++ ${PROBLEM}.cpp -y --open
submit/full:
	oj submit '${URL}' --language C++ ${PROBLEM}.cpp -y --open --full

score.txt: a.out tester.jar
	-rm score.txt
	for seed in $$(seq 1 100) ; do java -jar tester.jar -exec ./a.out -novis -seed $$seed | tee /dev/stderr | grep '{"seed":' >> score.txt ; done

score.txt/view:
	{ echo seed H W R C0 C P HW HW/R iteration time | tr ' ' '\t' ; cat score.txt | jq -r '"\(.seed)\t\(.H)\t\(.W)\t\(.R)\t\(.C0)\t\(.C)\t\(.P)\t\(.H*.W)\t\(.H*.W/.R|floor)\t\(.iteration)\t\(.time)"' ; } | sed 's/\t/\t| /g ; s/^/| / ; s/$$/\t|/' | expand -t 12 | sed '1 { p ; s/[^|]/-/g }'

standings:
	oj get-standings '${URL}' --format=csv | sed 's/\w\w\.\w\w\.\w\w\w\w \w\w:\w\w:\w\w/& EST/ ; y/ /~/ ; :1 ; s/,,/,-,/ ; t1 ; s/$$/,|/' | column -t -s , | sed 's/\( \+\)|\?/\1| /g ; s/^/| / ; s/.$$// ; y/~/ / ; 1 { p ; s/[^|]/-/g }'
