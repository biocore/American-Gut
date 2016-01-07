for file in `/bin/ls ${1}`
do
    if ! [[ ${file} =~ \.ipynb$ ]]; then
        continue
    fi

    runipy ${1}/${file} last_notebook
    if [ $? -ne 0 ]; then
        echo "Failed!"
        break
    fi
done
