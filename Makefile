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

score: a.out tester.jar
	for seed in $$(seq 1 20) ; do java -jar tester.jar -exec ./a.out -novis -seed $$seed | tee /dev/stderr | grep '{"seed":' >> score.txt ; done
	echo seed H W R C k delta score | tr ' ' '\t'
	cat score.txt | jq -r '"\(.seed)\t\(.H)\t\(.W)\t\(.R)\t\(.C)\t\(.k)\t\(.delta)\t\(.score)"'
