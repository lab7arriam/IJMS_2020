for name in $(cat bt_assemblys.tsv | awk -F\\t '{print $6}'); 
do       
	link2=$(echo $name | awk -F'/' '{print $10"_genomic.gbff.gz"}' | tr -d '\r');
		echo $name/$link2;
		wget $name/$link2;
done
