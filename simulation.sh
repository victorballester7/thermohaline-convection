#!/bin/bash

# Define the filename
filename="config/input.txt"

sumR_rho=0.02
MultiplyLe=10
IterR_rho=5
IterLe=5

startLe=0.01
startR_rho=0.96

# MultiplyR_rho=1
# MultiplyLe=2
# IterR_rho=1
# IterLe=10

# startLe=0.01
# startR_rho=1

editLe(){
    # Extract the value of Le
    Le=$(grep "^Le" "$filename" | awk '{print $2}')
    # Multiply the values by 10 (even in floating point arithmetic)
    new_Le=$(awk "BEGIN {print $Le * $MultiplyLe}")
    # Replace the lines with the new values
    sed -i "s/^Le $Le$/Le $new_Le/" "$filename"

}

editR_rho(){
    # Extract the value of R_rho
    R_rho=$(grep "^R_rho" "$filename" | awk '{print $2}')
    # Multiply the values by 10 (even in floating point arithmetic)
    new_R_rho=$(awk "BEGIN {print $R_rho + $sumR_rho}")
    # Replace the lines with the new values
    sed -i "s/^R_rho $R_rho$/R_rho $new_R_rho/" "$filename"
}

resetLe(){
    # Extract the value of Le
    Le=$(grep "^Le" "$filename" | awk '{print $2}')
    # Replace the lines with the new values
    sed -i "s/^Le $Le$/Le $startLe/" "$filename"
}

resetR_rho(){
    # Extract the value of R_rho
    R_rho=$(grep "^R_rho" "$filename" | awk '{print $2}')
    # Replace the lines with the new values
    sed -i "s/^R_rho $R_rho$/R_rho $startR_rho/" "$filename"
}

# ask if the user wants to remove the output folder
echo "Do you want to remove the output folder? (y/n)"
read answer
if [ "$answer" == "y" ]; then
    rm -rf output/Le_Rrho.txt
fi
resetLe
resetR_rho
for ((i=1; i<=$IterR_rho; i++)); do
    for ((j=1; j<=$IterLe; j++)); do
        echo $i/$IterR_rho $j/$IterLe
        ./run.sh 1
        editLe
    done
    resetLe
    editR_rho
done
resetR_rho