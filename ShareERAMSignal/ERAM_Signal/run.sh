set -e # Exit on error

peak_time="412"

while getopts ":t:" opt; do
    case $opt in
        t)
            peak_time="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            ;;
    esac
done

cd build
cmake ..
make -j 12
cd ..
time ./build/Signal $peak_time