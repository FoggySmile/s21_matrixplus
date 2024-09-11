CXX = g++
OS1LVL = $(shell uname)
OS2LVL = $(shell uname -o)
CXXFLAGS = -Wall -Werror -Wextra -g -std=c++17 -I/usr/local/include

ALL_FUNC = $(wildcard ./*.cpp)
SOURCE_TESTS_FUNC = $(wildcard ./tests/*.cpp ./tests/s21_matrix/*.cpp)

ifeq ($(OS1LVL), Darwin)
	det_OS = -lgtest -lgtest_main -lpthread
	GCOV_FLAGS = -fprofile-arcs -ftest-coverage
else
	ifeq ($(OS2LVL), GNU/Linux)
		det_OS = -lgtest -lgtest_main -lpthread -lrt -lsubunit
	else
		det_OS = -lgtest -lgtest_main -lpthread -lrt
		GCOV_FLAGS = -fprofile-arcs -ftest-coverage
	endif
endif

all: clean s21_matrix_oop.a test

install_gtest:
	sh gtest_installer.sh

test: 
	@$(CXX) $(CXXFLAGS) $(ALL_FUNC) $(SOURCE_TESTS_FUNC) $(GCOV_FLAGS) -g --coverage s21_matrix_oop.a -o s21_test_matrix $(det_OS)
	@./s21_test_matrix

s21_matrix_oop.a:
	@$(CC) $(CXXFLAGS) $(ALL_FUNC) -g -c
	@ar rc s21_matrix_oop.a *.o
	@ranlib s21_matrix_oop.a

gcov_report: clean s21_matrix_oop.a test
	@lcov -t "s21_test_matrix" -o test.info --no-external -c -d .
	@genhtml -o report test.info
	open report/index.html
	rm -rf *.gcno *gcda *.gco

clean:
	@rm -rf test *.gcda *.gcno *.o *.info *.a *.h.gch tests/*.gcno report 
	@rm -rf s21_test_matrix
	@rm -rf tests/*.h.gch 

style:
	@cp ../materials/linters/.clang-format .
	clang-format -n *.h $(ALL_FUNC) $(SOURCE_TESTS_FUNC)
	@rm -rf .clang-format

style_repair:
	@cp ../materials/linters/.clang-format .
	clang-format -i *.h $(ALL_FUNC) $(SOURCE_TESTS_FUNC)
	@rm -rf .clang-format

leaks:
	valgrind --tool=memcheck --leak-check=yes ./s21_test_matrix

install_brew:
	cd ~
	curl -fsSL https://rawgit.com/kube/42homebrew/master/install.sh | zsh
	brew install lcov
