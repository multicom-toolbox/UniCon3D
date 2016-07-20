

$num = @ARGV;
if($num != 2)
{
	die "The number of parameter is not correct! Para: <input> <Output>\n";
}

$input = $ARGV[0];
$output = $ARGV[1];
open(IN,"$input") || die "Failed to open file $input\n";
open(OUT,">$output") || die "Failed to write file $output\n";

while(<IN>)
{
	$line =$_;
	chomp $line;
	@content = split(/\s+/,$line);
	if(@content !=3)
	{
		next;
	}else{
		print OUT $content[0]." ".$content[1]." 0 8 ".$content[2]."\n";
	}
}
close IN;
close OUT;
