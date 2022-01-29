# How to build documentation

* install rust
* get the exhact code
	* The latest version of hte code can be found on the "main" branch.
	* either use git to clone a copy and switch to the relevant branch
	* or just download directly from github (you can choose the branch by clicking the dropdown menu towards the upper left of the screen)
* open a command shell and cd to the folder you downloaded
* run `cargo clean --doc` to clear any old documentation files
* build the docs
    * run `cargo doc --no-deps --open` to build the docs and automatically open in a browser window
    * run `cargo doc --no-deps` to just rebuild the docs (you'll have to refresh your browser to see the change)





