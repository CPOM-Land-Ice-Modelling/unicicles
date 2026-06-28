#!/bin/bash

CONFIG_DIR="configs"

for config_file in "$CONFIG_DIR"/*.yaml; do
    echo "Running with config: $config_file"
    python test_TerraFIRMA.py --config "$config_file" &
done
wait
echo "All tests completed."