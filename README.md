# De Bruijn Graph library

## Compilation

```bash
	mkdir -p build && cd build
	cmake .. && make
```

## Tests

* Run all tests
```bash
	mkdir -p build && cd build
	cmake -DCMAKE_BUILD_TYPE=DEBUG .. && make -j && ./tests/dbglib-tests
```

* Run a subset of tests
```bash
	mkdir -p build && cd build
	cmake -DCMAKE_BUILD_TYPE=DEBUG .. && make && ./tests/dbglib-tests --gtest_filter="regexp test names"
```   

For example, if you want to run test on the (virtual) Super-k-mer sorting, you can use
```bash
./tests/dbglib-tests --gtest_filter=*SkmerSorting*
```