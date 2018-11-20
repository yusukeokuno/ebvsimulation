#!/usr/local/bin/php
#$ -S /usr/local/bin/php
#$ -cwd
#$ -e log
#$ -o log
#$ -l s_vmem=4G,mem_req=4G
<?php

ini_set("memory_limit","2048M");
ini_set("error_reporting","E_ALL & ~E_WARNING & ~E_NOTICE");

//------------------------------------------------------------------------------
// simulation.php (for parallel calculation using SGE)
//------------------------------------------------------------------------------

function loadExons(){
	$in=fopen("genes.txt","r");
	$exons=array();
	while(!feof($in)){
		list($no,$posA,$posB,$type,$product_org,$product_edit)=explode("\t",rtrim(fgets($in)));
		if($type!=null){
			$exon=new StdClass;
			$exon->from=min($posA,$posB);
			$exon->to=max($posA,$posB);
			$exon->ID=$type.'@'.$product_edit;
			$exons[]=$exon;
		}
	}
	fclose($in);
	return $exons;
}

function getAffectedExons($from,$to,$exons){
	$affected=array();
	foreach($exons as $e){
		$add=false;
		// inside deletions
		if($from<=$e->from&&$e->to<=$to)$add=true;
		// edge of deletions
		if($e->from<=$from&&$from<=$e->to)$add=true;
		if($e->from<=$to&&$to<=$e->to)$add=true;
		//
		if($add)$affected[]=$e;
	}
	return $affected;
}

//------------------------------------------------------------------------------
// load deletions and calculate affected genes (observed number of deletions)
//------------------------------------------------------------------------------
$exons=loadExons();
$in=fopen("deletions.txt","r");
$out=fopen("out.txt","w");
$observed=array();
while(!feof($in)){
	$line=fgets($in);
	list($SampleID,$from,$to)=explode("\t",rtrim($line));
	if($to!=null){
		$affected=getAffectedExons($from,$to,$exons);
		$genes=array();
		foreach($affected as $a)$genes[$a->ID]=1;
		$line=rtrim($line)."\t";
		foreach($genes as $gene => $temp)$line.=$gene.",";
		$line.=PHP_EOL;
		fwrite($out,$line);
		foreach($genes as $gene => $temp)$observed[$gene]++;
	}
}
foreach($observed as $gene => $count)print $gene." ".$count.PHP_EOL;
fclose($out);

//------------------------------------------------------------------------------
// permutation
//------------------------------------------------------------------------------
$exons=loadExons();
$in=fopen("deletions.txt","r");
$lengths=array();
while(!feof($in)){
	$line=fgets($in);
	list($SampleID,$from,$to)=explode("\t",rtrim($line));
	if($to!=null)$lengths[]=$to-$from+1;
}

$temp=array();
foreach($exons as $exon)$temp[]=$exon->ID;
$allgenes=array_unique($temp);

// simulate for 1,000,000 times and compare with observed number of deletions
$out=fopen("result.".date("YmdHis").".".rand().".txt","w");
$score1=array();
$score2=array();
for($n=0;$n<1000000;$n++){
	$count=0;
	$simulated=array();
	foreach($lengths as $len){
		// place a simulated deletion
		$i=mt_rand(1,171323);
		$from=$i;
		$to=$i+$len-1;
		// get and count affected genes
		$a1=getAffectedExons($from,$to<171323?$to:171323,$exons);
		$a2=array();
		if($to>171323)$a2=getAffectedExons(1,$to-171323,$exons);
		$affected=array_merge($a1,$a2);
		$genes=array();
		foreach($affected as $a)$genes[$a->ID]=1;
		foreach($genes as $gene => $temp)$simulated[$gene]++;
	}
	// check if observed<=simulated (score1)
	foreach($allgenes as $gene)if($observed[$gene]<=$simulated[$gene])$score1[$gene]++;
	// check if observed>=simulated (score2)
	foreach($allgenes as $gene)if($observed[$gene]>=$simulated[$gene])$score2[$gene]++;
	//
	if($n%100==0)print $n." ";
}
fwrite($out,"gene\tobserved\tscore1\tscore2".PHP_EOL);
foreach($allgenes as $gene)fwrite($out,$gene."\t".(int)$observed[$gene]."\t".(int)$score1[$gene]."\t".(int)$score2[$gene].PHP_EOL);
fclose($out);

?>
