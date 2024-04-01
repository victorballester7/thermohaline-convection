#!/bin/bash

# Define the filename
filename="config/input.txt"

MultiplyR_rho=10
MultiplyLe=10
IterR_rho=6
IterLe=4

startLe=0.01
startR_rho=0.01

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
    new_R_rho=$(awk "BEGIN {print $R_rho * $MultiplyR_rho}")
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

for ((i=1; i<=$IterR_rho; i++)); do
    for ((j=1; j<=$IterLe; j++)); do
        ./run.sh 1
        editLe
    done
    resetLe
    editR_rho
done
resetR_rho