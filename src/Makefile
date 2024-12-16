FULL_FLAGS=-Werror -Wextra -Wall -std=c++11
ALL_SOURCES=$(shell find . -regextype posix-extended -regex ".*.(cc|h)")

all: s21_matrix_oop.a style cppcheck leaks gcov_report

s21_matrix_oop.a:
	g++ -c s21_matrix_oop.cc -o s21_matrix_oop.o
	ar rcs s21_matrix_oop.a s21_matrix_oop.o
	ranlib s21_matrix_oop.a

test:
	g++ s21_matrix_oop.cc cpp_test.cc -lgtest -fprofile-arcs -ftest-coverage -lgcov -o s21_test
	./s21_test

gcov_report: test
	lcov --no-external -t "s21_matrix" -o s21_matrix_oop.info -c -d . --ignore-errors mismatch 
	genhtml --ignore-errors mismatch -o report s21_matrix_oop.info
	open report/index.html

leaks: test
	leaks -atExit -- ./s21_test

valgrind: test
	valgrind ./s21_test

cppcheck:
	cppcheck --language=c++ --enable=all --suppress=missingIncludeSystem *.cc *.h

style:
	clang-format -n -style=google $(ALL_SOURCES)

style_fix:
	clang-format -i -style=google $(ALL_SOURCES)

add: style_fix
	git add $(ALL_SOURCES) Makefile

clean:
	-rm *.o *.a *.gcno *.gcda *.info s21_test
	-rm -rf report/
