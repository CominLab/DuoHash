.PHONY: build clean

build: build/libDuoHash.a

build/libDuoHash.a:
	mkdir -p build
	cd ./build && cmake .. && make

clean:
	rm -rf build
