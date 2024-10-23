#sudo apt install sra-toolkit 

mkdir -p data
list_id=(GSM4145661 GSM4145662 GSM4145663 GSM4145664 GSM4145665 GSM4145666)

for id in ${list_id[@]}
do
    fasterq-dump --threads 4 --progress $id -O data/
done

