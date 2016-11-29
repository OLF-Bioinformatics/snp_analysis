#!/usr/bin/env python

# python3

import sys
import csv
import xlsxwriter

filename = sys.argv[1].replace(".txt",".xlsx")
group = sys.argv[2]

wb = xlsxwriter.Workbook(filename)
ws = wb.add_worksheet(group)

with open(sys.argv[1],'r') as csvfile:
    table = csv.reader(csvfile, delimiter='\t')
    i = 0
    for row in table:
        ws.write_row(i, 0, row)
        i += 1

#Number of coplumns
col = len(row)


# Bold format
bold = wb.add_format({'bold': True})

# Background colors
formatA = wb.add_format({'bg_color':'#58FA82'})  # green
formatG = wb.add_format({'bg_color':'#F7FE2E'})  # yellow
formatC = wb.add_format({'bg_color':'#0000FF'})  # blue
formatT = wb.add_format({'bg_color':'#FF0000'})  # red
formatN = wb.add_format({'bg_color':'#E2CFDD'})  # pink
formatnormal = wb.add_format({'bg_color':'#FDFEFE'})  # white

# Font colors
formatlowqual = wb.add_format({'font_color':'#C70039', 'bg_color':'#E2CFDD'})  # red font on pink bg
formathighqual = wb.add_format({'font_color':'#000000', 'bg_color':'#FDFEFE'})  # blue font on white bg
formatperfectqual = wb.add_format({'font_color':'#0A028C', 'bg_color':'#FDFEFE'})
formatambigous = wb.add_format({'font_color':'#C70039', 'bg_color':'#E2CFDD'})  # red font on pink bg

# in formating: 1,2,3,4,
# 1 (row) and 2 (column) are first cell
# 3 (row) and 4 (column) are last cell
# Both rows and columns are zero indexed!
# Example: to higlight last row -> i-1,1,i-1,col-1

# order of conditions is very important
# once a condition is written for cell it cannot be over written

# Excel reads the qual vaules as text
# therefore cannot use numerical(type:cell with <) criteria
ws.conditional_format(i-1,1,i-1,col-1, {'type':'text',
                      'criteria':'containing',
                      'value':60,
                      'format':formatperfectqual})
ws.conditional_format(i-1,1,i-1,col-1, {'type':'text',
                      'criteria':'containing',
                      'value':59,
                      'format':formathighqual})
ws.conditional_format(i-1,1,i-1,col-1, {'type':'text',
                      'criteria':'not containing',
                      'value':100,
                      'format':formatlowqual})

ws.conditional_format(2,1,i-2,col-1, {'type':'cell',
                      'criteria':'==',
                      'value':'B$2',
                      'format':formatnormal})
ws.conditional_format(2,1,i-2,col-1, {'type':'text',
                      'criteria':'containing',
                      'value':'A',
                      'format':formatA})
ws.conditional_format(2,1,i-2,col-1, {'type':'text',
                      'criteria':'containing',
                      'value':'G',
                      'format':formatG})
ws.conditional_format(2,1,i-2,col-1, {'type':'text',
                      'criteria':'containing',
                      'value':'C',
                      'format':formatC})
ws.conditional_format(2,1,i-2,col-1, {'type':'text',
                      'criteria':'containing',
                      'value':'T',
                      'format':formatT})
ws.conditional_format(2,1,i-2,col-1, {'type':'text',
                      'criteria':'containing',
                      'value':'S',
                      'format':formatambigous})
ws.conditional_format(2,1,i-2,col-1, {'type':'text',
                      'criteria':'containing',
                      'value':'Y',
                      'format':formatambigous})
ws.conditional_format(2,1,i-2,col-1, {'type':'text',
                      'criteria':'containing',
                      'value':'R',
                      'format':formatambigous})
ws.conditional_format(2,1,i-2,col-1, {'type':'text',
                      'criteria':'containing',
                      'value':'W',
                      'format':formatambigous})
ws.conditional_format(2,1,i-2,col-1, {'type':'text',
                      'criteria':'containing',
                      'value':'K',
                      'format':formatambigous})
ws.conditional_format(2,1,i-2,col-1, {'type':'text',
                      'criteria':'containing',
                      'value':'M',
                      'format':formatambigous})
ws.conditional_format(2,1,i-2,col-1, {'type':'text',
                      'criteria':'containing',
                      'value':'N',
                      'format':formatN})

ws.set_column(0, 0, 30)
ws.set_column(1, col-1, 2)

# Add a freeze pane
ws.freeze_panes(2, 1)

# Rotate the header line by 90 degrees
format_rotation = wb.add_format({'rotation':'90'})
ws.set_row(0, 180, format_rotation)

#Format "mapping quality" cell
formatannotation = wb.add_format({'font_color':'#0A028C'})
ws.set_row(i-1, 15, formatannotation)

# Add bold to the "reference_call" cell
ws.conditional_format(1, 0, 1, 0, {'type':'text',
                      'criteria':'containing',
                      'value':'reference_call',
                      'format':bold})
# ws.set_row(1, None, bold)

wb.close()
