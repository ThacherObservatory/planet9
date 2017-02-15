# Planet 9 Repo

All of Thacher Astronomy's tools for searching for planet 9 are in this repository.

## Tasks
- [ ] Write up documentation
	- [ ] `gridder.py`
	- [ ] `make_stacked_image.py`
	-  [ ] `make_test_ims.py`
	-  [x] `pandastarrs.py`
	-  [ ] `planet9_movement.py`
	-  [ ] `planet9_reduction.py`
	-  [ ] `planet9_sim.py`
	-  [ ] `seeker.py`
	-  [x] `trylegal.py`
-  [x] automate film generation
-  [ ] automate planet9 recovery
	-  [ ] integrate background subtraction

## Index
* [gridder.py](#gridder.py)
* [make_stacked_image.py](#make_stacked_image.py)
* [make_test_ims.py](#make_test_ims.py)
* [pandastarrs.py](#pandastarrs.py)
* [planet9_movement.py](#planet9_movement.py)
* [planet9_reduction.py](#planet9_reduction.py)
* ~~planet9_reduction_ko.py~~
* [planet9_sim.py](#planet9_sim.py)
* [seeker.py](#seeker.py)
* [trylegal.py](#trylegal.py)

-----

### [gridder.py](https://github.com/ThacherObservatory/planet9/gridder.py)
### [make_stacked_image.py](https://github.com/ThacherObservatory/planet9/make_stacked_image.py)
### [make_test_ims.py](https://github.com/ThacherObservatory/planet9/make_test_ims.py)
### [pandastarrs.py](https://github.com/ThacherObservatory/planet9/pandastarrs.py)
Apparatus for extracting information from tables of PanSTARRS data.

| Functions | Arguments | Returns | Notes
| -------- | -------- | -------- | -------- |
|`load_info`||PANDAS dataframe|read PanSTARRS data into memory (defaults to `./data/PanSTARRS.table`)|
|`info_col`|string `column`| PANDAS dataframe|extracts only a particular column and its data from the PanSTARRS datatable
|`info_len`||`len` of data table|used to calculate the total number of entries in the PanSTARRS datafile|
### [planet9_movement.py](https://github.com/ThacherObservatory/planet9/planet9_movement.py)
### [planet9_reduction.py](https://github.com/ThacherObservatory/planet9/planet9_reduction.py)
Nick Edwards:
> it takes p9 images and subtracts the background and hcongrid stuff


### [planet9_sim.py](https://github.com/ThacherObservatory/planet9/planet9_sim.py)
Nick Edwards:
> it simulates p9 fields and makes movies


### [seeker.py](https://github.com/ThacherObservatory/planet9/seeker.py)
Nick Edwards:
> `seeker.py` is for automatically finding planet 9 and is unfinished

### [trylegal.py](https://github.com/ThacherObservatory/planet9/trylegal.py)
Apparatus for extracting information from tables of TRILEGAL data.

| Functions | Arguments | Returns | Notes
| -------- | -------- | -------- | -------- |
|`load_info`||PANDAS dataframe|read TRILEGAL data into memory (defaults to `./data/TRILEGAL.table`)|
|`info_col`|string `column`| PANDAS dataframe|extracts only a particular column and its data from the TRILEGAL datatable
|`info_len`||`len` of data table|used to calculate the total number of entries in the TRILEGAL datafile|