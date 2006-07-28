#!/usr/local/bin/perl

use strict;
use Spreadsheet::WriteExcel;

# Create a new Excel workbook called perl.xls
my $workbook = Spreadsheet::WriteExcel->new("/home/avb/Desktop/example1.xls");

# Add some worksheets
my $worksheet = $workbook->add_worksheet("Example");

# Add a Format
my $format = $workbook->add_format();
$format->set_bold();
$format->set_size(15);
$format->set_color('blue');
$format->set_align('center');

# Set the width of the first column in Worksheet
$worksheet->set_column(0, 0, 30);

# Set Worksheet as the active worksheet
$worksheet->activate();

# The general syntax is write($row, $col, $token, $format)

# Write some formatted text
$worksheet->write(0, 0, "Hello Excel!", $format);

my ($col, $row) = 0;
$worksheet->write($row, $col, "Hi Excel!", $format);
$worksheet->write(1,    $col, "Hi Excel!");

$worksheet->write('A3', 1.2345);
$worksheet->write('A4', '=SIN(PI()/4)');

# Write some unformatted text
$worksheet->write(2, 0, "One");
$worksheet->write(3, 0, "Two");

# Write some unformatted numbers
$worksheet->write(4, 0, 3);
$worksheet->write(5, 0, 4.00001);

# Write a number formatted as a date
my $date = $workbook->add_format();
$date->set_num_format('mmmm d yyyy h:mm AM/PM');
$worksheet->write(7, 0, 36050.1875, $date);

for my $i (0 .. 9) {
    my $value = $i*$i;
    $worksheet->write($i, 5, "$value"); # Cells F1 to F10
}
my $boldformat = $workbook->add_format();
$boldformat->set_bold();
$worksheet->write('F11', '=SUM(F1:F10)', $boldformat);

my @array  = (1,1.2,2.3,3.4);

$worksheet->write_row(1, 6, \@array);

my @eec =  (
         ['maggie', 'milly', 'molly', 'may'  ],
         [13,       14,      15,      16     ],
         ['shell',  'star',  'crab',  'stone']
        );

$worksheet->write_col('G9', \@eec);
$worksheet->write_row('G4', \@eec);


$worksheet->write_formula(0, 0, '=$B$3 + B4'  );
$worksheet->write_formula(1, 0, '=SIN(PI()/4)');
$worksheet->write_formula(2, 0, '=SUM(B1:B5)' );
$worksheet->write_formula('A4', '=IF(A3>1,"Yes", "No")'   );
$worksheet->write_formula('A5', '=AVERAGE(1, 2, 3, 4)'    );
$worksheet->write_formula('A6', '=DATEVALUE("1-Jan-2001")');

my $formula = $worksheet->store_formula('=F1 * 3 + 50');

for my $row (0..11) {
    $worksheet->repeat_formula($row, 6, $formula, $boldformat, 'F1', 'F'.($row +1));
}

my $textwrapformat = $workbook->add_format();
$textwrapformat->set_text_wrap();
$worksheet->write('D5', "It's\na bum\nwrap", $textwrapformat);

my $rotatedformat = $workbook->add_format();
$rotatedformat->set_rotation(30);
$worksheet->write('D8', "This text is rotated", $rotatedformat);

$workbook->close();
1;
