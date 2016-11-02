# Install Dependencies

1. [Assembly Dependencies](https://github.com/guyleonard/single_cell_workflow/tree/master/install_dependencies#1-assembly-dependencies)
2. [Gene Prediction Dependencies](https://github.com/guyleonard/single_cell_workflow/tree/master/install_dependencies#2-gene-prediction-dependencies)

## 1. Assembly Dependencies
A bash script to install the software:

## 2. Gene Prediction Dependencies
An [Ansible]() playbook to install the software:

This requires the user to obtain a username and password for access to [Repbase](http://www.girinst.org/repbase/), it is stored in an ansible 'vault' file.
This file is also password protected, so the RepeatMasker install will not work for user of this repo, you will need to make your own vault, containing your own password
with this command:

    ansible-vault create repbase_password.yml

and add your password like so:
```yaml
    ---
    repbase_password: PASSWORD
```

Your username is in the repeatmasker.yaml taskbook
