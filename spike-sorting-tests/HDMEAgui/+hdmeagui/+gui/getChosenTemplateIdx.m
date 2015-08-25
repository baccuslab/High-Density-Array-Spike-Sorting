function tIdx = getChosenTemplateIdx(handles)
    ah = handles.TemplateListbox;
    tIdx = get(ah, 'Value');
end